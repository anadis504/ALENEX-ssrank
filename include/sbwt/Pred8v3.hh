#ifndef _PRED8_V3_H_
#define _PRED8_V3_H_

#include <stdio.h>
#include <stdlib.h>

#include <chrono>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <ratio>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace std::chrono;

class Pred8v3 {
 public:
  Pred8v3() {}
  Pred8v3(const vector<uint64_t>& data) {
    _n = data.size();
    _min = data[0];
    _u = data[data.size() - 1] - _min;
    _nblocks = _u / 256 + 1;

    _X.resize(_nblocks + 1);                 // +1 for a useful dummy at the end
    /* _X = new uint32_t[_nblocks + 1];  */  //+1 for a useful dummy at the end
    for (uint64_t i = 0; i < _nblocks; i++) _X[i] = 0;
    uint8_t* _C = new uint8_t[_nblocks];
    for (uint64_t i = 0; i < _nblocks; i++) _C[i] = 0;

    // NB: note that because we substract _min from everything, the 0th bucket
    // is non-empty

    _nActiveBuckets = 0;
    for (uint64_t i = 0; i < _n; i++) {
      uint64_t v = data[i] - _min;
      if (!_X[v >> 8]) _nActiveBuckets++;
      _X[v >> 8]++;
    }

    // _Y = new uint8_t[_n];
    _Y.resize(_n);
    uint64_t yi = 0;
    for (uint64_t i = 0; i < _n;) {
      uint64_t v = data[i] - _min;
      uint64_t bcount = _X[v >> 8];
      _X[v >> 8] =
          (yi << 1) | (bcount > 0);  // LSB==1 indicates bucket is non-empty
      if (bcount) {
        _C[v >> 8] = (uint8_t)(bcount - 1);
        for (uint64_t j = 0; j < bcount; j++) {
          v = data[i] - _min;
          _Y[yi++] = v & 255;
          i++;
        }
      }  // else{
         //   cerr << "Should never happen: "<<i<<'\n';
         //}
    }
    uint32_t lastNonEmptyX = _X[0] & 0xFFFFFFFE;
    uint32_t lastNonEmptyC = _C[0];
    for (uint64_t i = 1; i < _nblocks; i++) {
      if (!_X[i]) {
        _X[i] = (((lastNonEmptyX >> 1) + lastNonEmptyC)
                 << 1);  // we want to point to it's last element
      } else {
        lastNonEmptyX = _X[i] & 0xFFFFFFFE;
        lastNonEmptyC = _C[i];
      }
    }
    _X[_nblocks] = (((_X[_nblocks - 1] >> 1) + _C[_nblocks - 1]) << 1);
    cerr << "_X[_nblocks]: " << _X[_nblocks] << '\n';
    // When answering select0(sq) the table _Z below can be used to shorten the
    // scan of upper-level block headers. In particular, _Z[(sq - _min)/4096]
    // gives the block number to safely start scanning from. If the density of
    // 1s is 2%, then the expected scan length is about 8 block headers.
    int64_t nz = _nblocks / 16;
    _Z.resize(nz);
    for (uint32_t i = 0; i < nz; i++) {
      _Z[i] = 0xFFFFFFFF;
    }
    for (uint64_t bi = 1; bi < _nblocks; bi++) {
      int64_t blockStart = bi << 8;
      int64_t onesBeforeBlock = (_X[bi] >> 1) - (_X[bi] & 1);
      int64_t zerosBeforeBlock =
          blockStart -
          onesBeforeBlock;  //+ _min; <--- NB: min not included in 0 count
      if (_Z[zerosBeforeBlock / 4096] == 0xFFFFFFFF) {
        _Z[zerosBeforeBlock / 4096] = bi - 1;  // 4096 = 256*16
      }
    }
    cerr << "Pred8v3: _u _n _min _nblocks _nActiveBuckets: " << _u << ' ' << _n
         << ' ' << _min << ' ' << _nblocks << ' ' << _nActiveBuckets << '\n';
    cerr << "Pred8v3: sizeInBytes(): " << sizeInBytes() << '\n';
    delete[] _C;
  }

  //  p is the index of the predecessor in the set
  //  bool is 1 if the value at index p is equal to key, else 0
  pair<int64_t, bool> inline getPred(int64_t key) const {
    if (key < _min) {
      return {-1, false};
    }
    key = key - _min;
    uint32_t x = _X[key >> 8];
    if (x & 1) {
      // non-empty bucket
      uint64_t y =
          x >> 1;  // this is where, in _Y, the bucket elements are located
      uint64_t bcount =
          ((_X[1 + (key >> 8)] >> 1) + (1 - (_X[1 + (key >> 8)] & 1))) -
          y;  //_C[key>>8] + 1;
      // uint64_t bcount = _C[key>>8] + 1;
      // cerr << "bcount: "<<bcount<<'\n';
      // cerr << "_C[key>>8]+1: "<<(_C[key>>8] + 1)<<'\n';
      uint64_t k = key & 255;
      for (uint64_t j = 0; j < bcount; j++) {
        if (_Y[y + j] >= k) {
          return {y + j - (1 - (_Y[y + j] == k)), (_Y[y + j] == k)};
        }
      }
      // getting here means that we did not find anything in the bucket that was
      // >= k this means the last element of the bucket is the predecessor of k,
      // and it is not equal to k
      return {y + bcount - 1, false};
    }
    // the bucket that key belongs to is empty
    return {(x >> 1), false};
  }

  int64_t rank(int64_t pos) const {
    pair<int64_t, bool> r = getPred(pos);
    return (r.first + 1) - ((uint64_t)(r.second));
  }

  int64_t select0(int64_t sq) const {
    if (sq < _min) {
      return sq;
    }
    if (sq > (_u + _min)) {
      return sq + _n;
    }
    int64_t sqq = sq - _min;
    // int64_t bi = sqq >> 8;
    int64_t bi = _Z[sqq / 4096LL];  // 4096 = 256*16 = 2^12
    /* if (sqq >= 4235558912) {
      cerr << "select0: sq = " << sq << " sqq = " << sqq << " bi = " << bi
           << '\n';
    } */

    for (; bi < _nblocks; bi++) {
      int64_t blockStart = bi << 8;
      int64_t nextBlockStart = (bi + 1) << 8;
      int64_t onesBeforeBlock = _X[bi] >> 1;
      int64_t onesBeforeNextBlock = _X[1 + bi] >> 1;
      int64_t zerosBeforeBlock = blockStart - onesBeforeBlock + _min;
      int64_t zerosBeforeNextBlock =
          nextBlockStart - onesBeforeNextBlock + _min;
      if (zerosBeforeNextBlock >= sq) {
        // we've found our block
        int64_t pos = blockStart;
        if (_X[bi] & 1) {
          // block is not empty, we need to scan it
          int64_t y = _X[bi] >> 1;  // this is where the elements start in _Y
          uint64_t bcount = ((_X[1 + bi] >> 1) + (1 - (_X[1 + bi] & 1))) -
                            y;  // bcount is # elements in block bi
          int64_t prev = 0;
          int64_t zerosBeforePos = zerosBeforeBlock;
          for (uint64_t j = 0; j < bcount; j++) {
            int64_t v = _Y[y + j];
            int64_t zeroesBetween = v - prev;
            if (zerosBeforePos + zeroesBetween >= sq) {
              return pos + (sq - zerosBeforePos);
            }
            zerosBeforePos += zeroesBetween;
            prev = v + 1;
            pos = blockStart + v + 1;
          }
          // getting here means that the sq^th zero is between the last element
          // of this block and the start of the next block
          return pos + (sq - zerosBeforePos);
        }
        // getting here means the block that contains the sq^th 0 is empty
        // do some arithemtic and bail
        return blockStart + (sq - zerosBeforeBlock) + 1;
      }
      /* uint64_t up = onesBeforeBlock >> 8;
      bi += (onesBeforeBlock >> 8) ? (onesBeforeBlock >> 8) - 1
                                   : 0;  */ // 4096 = 256*16 = 2^12
    }
    // Should not happen I suppose
    return 0;
  }

  int64_t select0_slow(int64_t sq) const {
    if (sq < _min) {
      return sq;
    }
    if (sq > (_u + _min)) {
      return sq + _n;
    }
    /* if (sq > 410510) {
      cout << "Searching for sq = " << sq << '\n';
    } */
    int64_t sqq = sq - _min;
    int64_t bi = sqq >> 8;
    // int64_t bi = _Z[sqq / 4096];  // 4096 = 256*16 = 2^12
    cout << "select0_slow: sq = " << sq << " sqq = " << sqq << " bi = " << bi
         << '\n';
    for (; bi < _nblocks; bi++) {
      int64_t blockStart = bi << 8;
      int64_t nextBlockStart = (bi + 1) << 8;
      int64_t onesBeforeBlock = _X[bi] >> 1;
      int64_t onesBeforeNextBlock = _X[1 + bi] >> 1;
      int64_t zerosBeforeBlock = blockStart - onesBeforeBlock + _min;
      int64_t zerosBeforeNextBlock =
          nextBlockStart - onesBeforeNextBlock + _min;
      if (zerosBeforeNextBlock >= sq) {
        cout << "select0_slow: found block bi = " << bi
             << " blockStart = " << blockStart
             << " nextBlockStart = " << nextBlockStart
             << " zerosBeforeBlock = " << zerosBeforeBlock
             << " zerosBeforeNextBlock = " << zerosBeforeNextBlock << '\n';
        // we've found our block
        int64_t pos = blockStart;
        if (_X[bi] & 1) {
          // block is not empty, we need to scan it
          int64_t y = _X[bi] >> 1;  // this is where the elements start in _Y
          uint64_t bcount = ((_X[1 + bi] >> 1) + (1 - (_X[1 + bi] & 1))) -
                            y;  // bcount is # elements in block bi
          int64_t prev = 0;
          int64_t zerosBeforePos = zerosBeforeBlock;
          for (uint64_t j = 0; j < bcount; j++) {
            int64_t v = _Y[y + j];
            int64_t zeroesBetween = v - prev;
            if (zerosBeforePos + zeroesBetween >= sq) {
              return pos + (sq - zerosBeforePos);
            }
            zerosBeforePos += zeroesBetween;
            prev = v + 1;
            pos = blockStart + v + 1;
          }
          // getting here means that the sq^th zero is between the last element
          // of this block and the start of the next block
          return pos + (sq - zerosBeforePos);
        }
        // getting here means the block that contains the sq^th 0 is empty
        // do some arithemtic and bail
        return blockStart + (sq - zerosBeforeBlock) + 1;
      }
    }
    // Should not happen I suppose
    return 0;
  }

  size_t getu() const { return _u; }
  size_t getn() const { return _n; }

  uint64_t sizeInBytes() const {
    return sizeof(_u) + sizeof(_n) + sizeof(_min) + sizeof(_nblocks) +
           _X.size() * sizeof(uint32_t) + _Y.size() +
           _Z.size() * sizeof(uint32_t);
  }

  int64_t serialize(std::ostream& os) const {
    uint64_t written = 0;

    os.write(reinterpret_cast<const char*>(&_u), sizeof(_u));
    os.write(reinterpret_cast<const char*>(&_n), sizeof(_n));
    os.write(reinterpret_cast<const char*>(&_min), sizeof(_min));
    os.write(reinterpret_cast<const char*>(&_nblocks), sizeof(_nblocks));

    uint64_t zsize = _Z.size();
    cout << "Pred8v3::serialize: _X.size() = "
         << (_nblocks + 1) * sizeof(uint32_t) << " _Y.size() = " << _n
         << " _Z.size() = " << zsize * sizeof(uint32_t) << '\n';

    os.write(reinterpret_cast<const char*>(_X.data()),
             (_nblocks + 1) * sizeof(uint32_t));

    os.write(reinterpret_cast<const char*>(_Y.data()), _n * sizeof(uint8_t));

    os.write(reinterpret_cast<const char*>(&zsize), sizeof(zsize));
    os.write(reinterpret_cast<const char*>(_Z.data()),
             zsize * sizeof(uint32_t));

    written += 4 * sizeof(uint64_t) + sizeof(uint64_t) +
               (_nblocks + 1) * sizeof(uint32_t) + _n +
               zsize * sizeof(uint32_t);

    return written;
  }

  void load(std::istream& is) {
    is.read(reinterpret_cast<char*>(&_u), sizeof(_u));
    is.read(reinterpret_cast<char*>(&_n), sizeof(_n));
    is.read(reinterpret_cast<char*>(&_min), sizeof(_min));
    is.read(reinterpret_cast<char*>(&_nblocks), sizeof(_nblocks));

    uint64_t zsize;

    _X.resize(_nblocks + 1);
    is.read(reinterpret_cast<char*>(_X.data()),
            (_nblocks + 1) * sizeof(uint32_t));

    _Y.resize(_n);
    is.read(reinterpret_cast<char*>(_Y.data()), _n * sizeof(uint8_t));

    is.read(reinterpret_cast<char*>(&zsize), sizeof(zsize));
    _Z.resize(zsize);
    is.read(reinterpret_cast<char*>(_Z.data()), zsize * sizeof(uint32_t));
  }

  ~Pred8v3() = default;

  Pred8v3(Pred8v3& other) {
    this->_u = other._u;
    this->_n = other._n;
    this->_min = other._min;
    this->_nblocks = other._nblocks;
    this->_nActiveBuckets = other._nActiveBuckets;
    this->_X = other._X;
    this->_Y = other._Y;
    this->_Z = other._Z;
  }

 private:
  uint64_t _u = 0;    // universe size
  uint64_t _n = 0;    // number of elements
  uint64_t _min = 0;  // value of the smallest element
  uint64_t _nblocks = 0;
  uint64_t _nActiveBuckets = 0;
  std::vector<uint32_t> _X;
  std::vector<uint8_t> _Y;
  std::vector<uint32_t> _Z;
};

#endif
