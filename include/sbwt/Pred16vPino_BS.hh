#ifndef _PRED16_PINO_BS_H_
#define _PRED16_PINO_BS_H_

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

#include "Pred8vS1_new.hh"

using namespace std;
using namespace std::chrono;

class Pred16Pino_BS {
 public:
  Pred16Pino_BS() {}
  Pred16Pino_BS(const vector<uint64_t>& data) {
    _n = data.size();
    _min = data[0];
    _u = data[data.size() - 1] - _min;
    _nblocks = (_u >> 16) + ((_u % 65536) > 1);

    _X.resize(_nblocks + 1);
    // _X = new uint32_t[_nblocks + 1];  //+1 for a useful dummy at the end
    for (uint64_t i = 0; i < _nblocks; i++) _X[i] = 0;
    uint16_t* _C = new uint16_t[_nblocks];
    for (uint64_t i = 0; i < _nblocks; i++) _C[i] = 0;

    // NB: note that because we substract _min from everything, the 0th bucket
    // is non-empty

    _nActiveBuckets = 0;
    for (uint64_t i = 0; i < _n; i++) {
      uint64_t v = data[i] - _min;
      if (!_X[v >> 16]) _nActiveBuckets++;
      _X[v >> 16]++;
    }

    // _Y = new uint16_t[_n];
    _Y.resize(_n);

    cerr << "Pred16Pino_BS: _u _n _min _nblocks _nActiveBuckets: " << _u << ' '
         << _n << ' ' << _min << ' ' << _nblocks << ' ' << _nActiveBuckets
         << " _nblocks: " << _nblocks << " _nblocks >> 16 " << (_nblocks << 16)
         << " _X[_nblocks] " << _X[_nblocks] << '\n';
    uint64_t yi = 0;
    for (uint64_t i = 0; i < _n;) {
      uint64_t v = data[i] - _min;
      uint64_t bcount = _X[v >> 16];
      _X[v >> 16] = yi;
      //(yi << 1) | (bcount > 0);  // LSB==1 indicates if bucket is non-empty
      if (bcount) {
        _C[v >> 16] = (uint16_t)(bcount - 1);
        for (uint64_t j = 0; j < bcount; j++) {
          v = data[i] - _min;
          _Y[yi++] = v & 65535;
          i++;
        }
      }  // else{
         //   cerr << "Should never happen: "<<i<<'\n';
         //}
    }
    cerr << "_X[_nblocks]: " << _X[_nblocks] << '\n';
    
    if (!_X[_nblocks]) _X[_nblocks] = _X[_nblocks-1];
    cerr << "_X[_nblocks]: " << _X[_nblocks] << '\n';
    uint32_t prevNonEmptyX = _X[_nblocks];
    for (uint64_t i = _nblocks; i > 0; i--) {
      if (!_X[i]) {
        _X[i] = prevNonEmptyX;  // we want to point to it's last element
      } else {
        prevNonEmptyX = _X[i];
      }
    }
    //_X[_nblocks] = (((_X[_nblocks - 1] >> 1) + _C[_nblocks - 1]) << 1);

    delete[] _C;

    std::vector<uint64_t> block_indices;
    // the escape index as the Pred8 structure always starts from an element
    uint64_t bi = 0;
    block_indices.push_back(bi);
    for (uint64_t i = 0; i < _nblocks; i++) {
      uint32_t x = _X[i];
      uint64_t y = x;
      uint64_t bcount = (_X[1 + i]) ? (_X[1 + i]) - y : 0;
      bi += 1 + bcount;
      block_indices.push_back(bi);
    }
    _U = Pred8vS1_new(block_indices);

    cout << "returned from Pred8vS1_new constructor" << endl;
    uint64_t wrongs = 0;
    /* for (uint64_t bi = 0; bi < _nblocks; bi++) {
      uint32_t x = _X[bi];
      uint64_t y = x;
      uint64_t bcount = _X[1 + bi] - y;
      auto _y = _U.select1(bi);

      // the returned element should point to the first position of the block
      // or, if block is empty to its predecessor
      if (y != _y.first) {
        wrongs++;
        cout << "Pred8vPinoLight: ERROR: bi = " << bi << " y = " << y
        << " _y.first = " << _y.first << " block empty = " << !(x & 1)
        << '\n';
        }
        if (bcount != _y.second) {
          wrongs++;
          cout << "Pred8vPinoLight: ERROR: bcount != _y.second bi = " << bi
          << " bcount = " << bcount << " _y.second = " << _y.second
          << " block empty = " << !bcount << '\n';
          }

          if (wrongs > 50) {
            exit(1);
            }
            } */
    cerr << "Pred16Pino_BS: sizeInBytes(): " << sizeInBytes() << '\n';
  }

  //  p is the index of the predecessor in the set
  //  bool is 1 if the value at index p is equal to key, else 0
  pair<int64_t, bool> inline getPred(int64_t key) const {
    if (key < _min) {
      return {-1, false};
    }
    key = key - _min;
    uint32_t x = _X[key >> 16];
    if (x & 1) {
      // non-empty bucket
      uint64_t y =
          x >> 1;  // this is where, in _Y, the bucket elements are located,
                   // starting with the # of items in the bucket
      uint64_t bcount =
          ((_X[1 + (key >> 16)] >> 1) + (1 - (_X[1 + (key >> 16)] & 1))) -
          y;  //_C[key>>8] + 1;
      // uint64_t bcount = _C[key>>8] + 1;
      // cerr << "bcount: "<<bcount<<'\n';
      // cerr << "_C[key>>8]+1: "<<(_C[key>>8] + 1)<<'\n';
      uint64_t k = key & 65535;
      if (bcount > (1 << 6)) {
        uint64_t count = bcount;
        uint64_t it, step, j;
        j = y;
        while (count > 0) {
          it = j;
          step = count / 2;
          it += step;

          if (_Y[it] < k) {
            j = ++it;
            count -= step + 1;
          } else
            count = step;
        }
        return {j - (1 - (_Y[j] == k)), (_Y[j] == k)};
      } else {
        for (uint64_t j = 0; j < bcount; j++) {
          if (_Y[y + j] >= k) {
            return {y + j - (1 - (_Y[y + j] == k)), (_Y[y + j] == k)};
          }
        }
        // getting here means that we did not find anything in the bucket that
        // was >= k this means the last element of the bucket is the predecessor
        // of k, and it is not equal to k
        return {y + bcount - 1, false};
      }
    }
    // the bucket that key belongs to is empty
    return {(x >> 1), false};
  }

  pair<int64_t, bool> getPred_Select(int64_t key) const {
    if (key < _min) {
      return {-1, false};
    }
    key = key - _min;
    int64_t bi = key >> 16;
    auto p = _U.select1(bi);

    uint64_t y = p.first;
    int64_t bcount = p.second;
    // cout << " Called rank for key " << key << " bi " << bi << " first block
    // id "
    //      << y << " bcount " << bcount << endl;
    if (bcount) {
      uint64_t k = key & 65535;
      if (bcount > (1 << 6)) {
        uint64_t count = bcount;
        uint64_t it, step, j;
        j = y;
        while (count > 0) {
          it = j;
          step = count / 2;
          it += step;

          if (_Y[it] < k) {
            j = ++it;
            count -= step + 1;
          } else
            count = step;
        }
        return {j - (1 - (_Y[j] == k)), (_Y[j] == k)};
      } else {
        for (uint64_t j = 0; j < bcount; j++) {
          // cout << " _Y[y + j] " << _Y[y + j] << endl;
          if (_Y[y + j] >= k) {
            // cout << " returning y + j - (1 - (_Y[y + j] == k)), (_Y[y + j] ==
            // k) "
            // << (y + j - (1 - (_Y[y + j] == k))) << endl;
            return {y + j - (1 - (_Y[y + j] == k)), (_Y[y + j] == k)};
          }
        }
      }
      // getting here means that we did not find anything in the bucket that
      // was >= k this means the last element of the bucket is the predecessor
      // of k, and it is not equal to k
      return {y + bcount - 1, false};
    }
    return {y, false};
  }

  /* int64_t rank(int64_t pos) const {
    pair<int64_t, bool> r = getPred(pos);
    return (r.first + 1) - ((uint64_t)(r.second));
  } */

  int64_t rank(int64_t pos) const {
    pair<int64_t, bool> r = getPred_Select(pos);
    return (r.first + 1) - ((uint64_t)(r.second));
  }

  size_t getu() const { return _u; }
  size_t getn() const { return _n; }

  uint64_t sizeInBytes() const {
    uint64_t sz =
        4 * sizeof(uint64_t) + _U.sizeInBytes() + (sizeof(uint16_t) * _n);
    return sz;
  }

  int64_t serialize(std::ostream& os) const {
    uint64_t written = 0;

    os.write(reinterpret_cast<const char*>(&_u), sizeof(_u));
    os.write(reinterpret_cast<const char*>(&_n), sizeof(_n));
    os.write(reinterpret_cast<const char*>(&_min), sizeof(_min));
    os.write(reinterpret_cast<const char*>(&_nblocks), sizeof(_nblocks));

    cout << "Pred16vPino_BS::serialize: _Y.size() = " << _n
         << " _U sizeInBytes " << _U.sizeInBytes() << '\n';

    /* os.write(reinterpret_cast<const char*>(_X.data()),
             (_nblocks + 1) * sizeof(uint32_t)); */

    os.write(reinterpret_cast<const char*>(_Y.data()), _n * sizeof(uint8_t));

    written += 4 * sizeof(uint64_t) + sizeof(uint64_t) + _n;
    written += _U.serialize(os);
    return written;
  }
  void load(std::istream& is) {
    is.read(reinterpret_cast<char*>(&_u), sizeof(_u));
    is.read(reinterpret_cast<char*>(&_n), sizeof(_n));
    is.read(reinterpret_cast<char*>(&_min), sizeof(_min));
    is.read(reinterpret_cast<char*>(&_nblocks), sizeof(_nblocks));

    /* _X.resize(_nblocks + 1);
    is.read(reinterpret_cast<char*>(_X.data()),
            (_nblocks + 1) * sizeof(uint32_t)); */

    _Y.resize(_n);
    is.read(reinterpret_cast<char*>(_Y.data()), _n * sizeof(uint8_t));
    ;
    _U.load(is);
  }

  Pred16Pino_BS(Pred16Pino_BS& other) {
    this->_u = other._u;
    this->_n = other._n;
    this->_min = other._min;
    this->_nblocks = other._nblocks;
    this->_nActiveBuckets = other._nActiveBuckets;
    this->_X = other._X;
    this->_Y = other._Y;
    this->_U = other._U;
  }

 private:
  uint64_t _u = 0;    // universe size
  uint64_t _n = 0;    // number of elements
  uint64_t _min = 0;  // value of the smallest element
  uint64_t _nblocks = 0;
  uint64_t _nActiveBuckets = 0;
  /* uint32_t* _X;
  uint16_t* _Y; */
  // uint8_t *_C;
  std::vector<uint32_t> _X;
  std::vector<uint16_t> _Y;
  Pred8vS1_new _U;
};

#endif
