#ifndef _PRED8_VS1_H_
#define _PRED8_VS1_H_

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

class Pred8vS1 {
 public:
  Pred8vS1() {}
  Pred8vS1(const vector<uint64_t>& data) {
    cout << "Start of Pred8vS1 constructor\n";
    _n = data.size();
    _min = data[0];
    _u = data[data.size() - 1] - _min;
    _nblocks = _u / 256 + ((_u % 256) > 1);

    _X.resize(_nblocks + 1);  //+1 for a useful dummy at the end
    for (uint64_t i = 0; i < _nblocks; i++) _X[i] = 0;
    uint8_t* _C = new uint8_t[_nblocks];
    for (uint64_t i = 0; i < _nblocks; i++) _C[i] = 0;

    // NB: note that because we substract _min from everything, the 0th bucket
    // is non-empty
    // NEW NB: note that because we are compressing the top level layer of a
    // Predv8, the first element there will be a 0. Technically, min is not
    // needed

    _nActiveBuckets = 0;
    for (uint64_t i = 0; i < _n; i++) {
      uint64_t v = data[i] - _min;
      if (!_X[v >> 8]) _nActiveBuckets++;
      _X[v >> 8]++;
    }
    uint64_t n_hints = _n / 32 + 1;
    _Y.resize(_n + n_hints * 4);  // +4 bytes per hint

    cerr << "Pred8vS1: _u _n _min _nblocks _nActiveBuckets: " << _u << ' ' << _n
         << ' ' << _min << ' ' << _nblocks << ' ' << _nActiveBuckets << '\n';
    cerr << "Pred8vS1: sizeInBytes(): " << sizeInBytes() << '\n';
    uint64_t yi = 0;
    uint64_t y_ptr = 0;
    for (uint64_t i = 0; i < _n;) {
      uint64_t v = data[i] - _min;
      uint64_t bcount = _X[v >> 8];
      _X[v >> 8] =
          (yi << 1) | (bcount > 0);  // LSB==1 indicates if bucket is non-empty
      if (bcount) {
        _C[v >> 8] = (uint8_t)(bcount - 1);
        for (uint64_t j = 0; j < bcount; j++) {
          if (yi % 32 == 0) {
            ((uint32_t*)(_Y.data() + y_ptr))[0] = (uint32_t)(v >> 8);
            y_ptr += 4;
          }
          v = data[i] - _min;
          _Y[y_ptr] = v & 255;
          y_ptr++;
          yi++;
          i++;
        }
      }  // else{
      //   cerr << "Should never happen: "<<i<<'\n';
      //}
    }
    if (_Y.size() != y_ptr) {
      cout << " _Y.size() " << _Y.size() << " y_ptr " << y_ptr << endl;
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
    cerr << "(_X[_nblocks]>>1): " << (_X[_nblocks] >> 1) << '\n';
    delete[] _C;
  }

  // This select assumes an unoverse of _u/256 + m with no repeated elements
  // it returns what the _X[sq << 8] would return
  // and the number of elements in the bucket
  pair<int64_t, uint32_t> select1(int64_t sq) const {
    // the position in _Y of the first byte of our relevant hint
    uint64_t h_idx = (sq / 32) * 36;
    uint64_t el_idx = h_idx + 4 + (sq % 32);
    uint64_t el = _Y[el_idx];
    uint64_t next_el_idx = ((sq + 1) / 32) * 36 + 4 + ((sq + 1) % 32);
    uint64_t next_el = _Y[next_el_idx];
    // uint32_t hint = ((uint32_t*)((_Y.data() + h_idx)))[0];
    uint32_t hint = *reinterpret_cast<const uint32_t*>(_Y.data() + h_idx);
    /* if (sq >= 1363663)
      cout << "Hint for " << sq << " at " << h_idx << " : " << hint
           << " el_ind " << el_idx << endl; */
    uint32_t x = _X[hint];
    // non-empty bucket
    uint64_t y = x >> 1;  // this is where, in _Y, the bucket elements are
                          // located, starting with the # of items in the bucket
    uint64_t j = 1;
    uint64_t next_y = _X[hint + j] >> 1;

    // we need to know when to stop scanning block so next_y has to point to
    // strictly after sq
    while (next_y <= sq && (hint + j) < _nblocks) {
      j++;
      next_y = _X[hint + j] >> 1;
    }
    // now we scan from the upper bucket of our current element to find the
    // upper bucket of the next element
    uint64_t j_2 = j - 1;  // this is the bucket of the curr el
    uint64_t next_y_2 = _X[hint + j_2] >> 1;
    while (next_y_2 <= (sq + 1) && (hint + j_2) < _nblocks) {
      j_2++;
      next_y_2 = _X[hint + j_2] >> 1;
    }
    uint64_t upperbits = (hint + j - 1) << 8;
    uint64_t next_upperbits = ((hint + j_2 - 1) << 8);

    uint64_t curr_pos = upperbits + (uint64_t)el;
    uint64_t next_pos = next_upperbits + (uint64_t)next_el;
    // bcount will never get below zero as there are no repeated elements in our
    // set
    uint64_t bcount = (next_pos - curr_pos) - 1;

    // update position to be returned as we don't want to include the count of
    // 0's.
    curr_pos -= sq;
    // if the block is empty, the position must point to the predecessor of the
    // block
    curr_pos -= (bool)(bcount == 0);

    /* if (sq == 2521747 or sq == 450198 or sq == 653898) {
      cout << "sq " << sq << " hint + j - 1: " << hint + j - 1 << " upperbits "
           << upperbits << " el " << el << " next upperbits " << next_upperbits
           << " result - sq: " << curr_pos << " next_el " << next_el << endl;
    } */

    return {(int64_t)curr_pos, (uint32_t)bcount};
  }

  int64_t select1_scanning(int64_t sq) const {
    // the position in _Y of the first byte of our relevant hint
    uint64_t h_idx = (sq / 32) * 36;
    uint64_t el_idx = h_idx + 4 + (sq % 32);
    uint64_t el = _Y[el_idx];
    // uint32_t hint = ((uint32_t*)((_Y.data() + h_idx)))[0];

    uint32_t x = _X[0];
    // non-empty bucket
    uint64_t y = x >> 1;  // this is where, in _Y, the bucket elements are
                          // located, starting with the # of items in the bucket
    uint64_t j = 1;
    uint64_t next_y = _X[j] >> 1;
    while (next_y < sq && j < _nblocks) {
      j++;
      uint64_t next_y = _X[j] >> 1;
      cout << "select1_scanning: j = " << j << " next_y = " << next_y << endl;
    }

    uint64_t upperbits = (j - 1) << 8;
    uint64_t result = upperbits + (uint64_t)el;
    return (result + _min);
  }

  size_t getu() const { return _u; }
  size_t getn() const { return _n; }

  uint64_t sizeInBytes() const {
    uint64_t sz = 4 * sizeof(uint64_t) + (sizeof(uint32_t) * (_nblocks + 1)) +
                  _n + _n / 32 * 4;
    return sz;
  }

  int64_t serialize(std::ostream& os) const {
    uint64_t written = 0;

    os.write(reinterpret_cast<const char*>(&_u), sizeof(_u));
    os.write(reinterpret_cast<const char*>(&_n), sizeof(_n));
    os.write(reinterpret_cast<const char*>(&_min), sizeof(_min));
    os.write(reinterpret_cast<const char*>(&_nblocks), sizeof(_nblocks));

    cout << "Pred8vS1::serialize: _X.size() = "
         << (_nblocks + 1) * sizeof(uint32_t) << " _Y.size() = " << _n << '\n';

    os.write(reinterpret_cast<const char*>(_X.data()),
             (_nblocks + 1) * sizeof(uint32_t));

    uint64_t n_hints = _n / 32 + 1;
    os.write(reinterpret_cast<const char*>(_Y.data()),
             (_n + n_hints * 4) * sizeof(uint8_t));

    written += 4 * sizeof(uint64_t) + sizeof(uint64_t) +
               _X.size() * sizeof(uint32_t) + _Y.size() * sizeof(uint8_t);

    return written;
  }

  void load(std::istream& is) {
    is.read(reinterpret_cast<char*>(&_u), sizeof(_u));
    is.read(reinterpret_cast<char*>(&_n), sizeof(_n));
    is.read(reinterpret_cast<char*>(&_min), sizeof(_min));
    is.read(reinterpret_cast<char*>(&_nblocks), sizeof(_nblocks));

    _X.resize(_nblocks + 1);
    is.read(reinterpret_cast<char*>(_X.data()),
            (_nblocks + 1) * sizeof(uint32_t));
    uint64_t n_hints = _n / 32 + 1;
    _Y.resize(_n + n_hints * 4);  // +4 bytes per hint
    is.read(reinterpret_cast<char*>(_Y.data()),
            (_n + n_hints * 4) * sizeof(uint8_t));
  }

  Pred8vS1(Pred8vS1& other) {
    this->_u = other._u;
    this->_n = other._n;
    this->_min = other._min;
    this->_nblocks = other._nblocks;
    this->_nActiveBuckets = other._nActiveBuckets;
    this->_X = other._X;
    this->_Y = other._Y;
  }

 private:
  uint64_t _u = 0;    // universe size
  uint64_t _n = 0;    // number of elements
  uint64_t _min = 0;  // value of the smallest element
  uint64_t _nblocks = 0;
  uint64_t _nActiveBuckets = 0;
  std::vector<uint32_t> _X;
  std::vector<uint8_t> _Y;
};

#endif
