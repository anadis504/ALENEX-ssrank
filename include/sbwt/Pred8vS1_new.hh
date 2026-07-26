#ifndef _PRED8_VS1_NEW_H_
#define _PRED8_VS1_NEW_H_

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

class Pred8vS1_new {
 public:
  Pred8vS1_new() {}
  Pred8vS1_new(const vector<uint64_t>& data) {
    //cout << "Start of Pred8vS1_new constructor\n";
    _n = data.size();
    _u = data[data.size() - 1];
    _nblocks = _u / 256 + ((bool)(_u % 256) > 0);

    _X.resize(_nblocks + 1);  //+1 for a useful dummy at the end
    for (uint64_t i = 0; i < _nblocks; i++) _X[i] = 0;

    // NB: note that because we substract _min from everything, the 0th bucket
    // is non-empty
    // NEW NB: note that because we are compressing the top level layer of a
    // Predv8, the first element there will be a 0. Technically, min is not
    // needed

    _nActiveBuckets = 0;
    for (uint64_t i = 0; i < _n; i++) {
      uint64_t v = data[i];
      assert(data[i] <= _u && data[i] >= 0);
      if (!_X[v >> 8]) _nActiveBuckets++;
      _X[v >> 8]++;
    }
    uint64_t n_hints = _n / 32 + 1;
    /* cerr << "Pred8vS1_new: _u _n _nblocks _nActiveBuckets: " << _u << ' '
         << _n << ' ' << ' ' << _nblocks << ' ' << _nActiveBuckets
         << '\n'; */
    _Y.resize(_n + n_hints * 4);  // +4 bytes per hint

    /* cerr << " _nblocks: " << _nblocks << " n_hints: " << n_hints
         << " _nblocks << 8 " << (_nblocks << 8) << " _X[_nblocks] "
         << _X[_nblocks] << '\n';
    cerr << "Pred8vS1_new: sizeInBytes(): " << sizeInBytes() << '\n'; */
    uint64_t yi = 0;
    uint64_t y_ptr = 0;
    for (uint64_t i = 0; i < _n;) {
      uint64_t v = data[i];
      uint64_t bcount = _X[v >> 8];
      _X[v >> 8] = yi;
      /* (yi << 1) | (bcount > 0) */
      if (bcount) {
        for (uint64_t j = 0; j < bcount; j++) {
          if (yi % 32 == 0) {
            ((uint32_t*)(_Y.data() + y_ptr))[0] = (uint32_t)(v >> 8);
            y_ptr += 4;
          }
          v = data[i];
          _Y[y_ptr] = v & 255;
          y_ptr++;
          yi++;
          i++;
        }
      }  
    }
    if (_Y.size() != y_ptr) {
      cout << " _Y.size() " << _Y.size() << " y_ptr " << y_ptr << endl;
    }
    
    uint32_t prevNonEmptyX = _X[_nblocks - 1];
    // uint32_t lastNonEmptyC = _C[0];
    for (uint64_t i = _nblocks - 1; i > 0; i--) {
      if (!_X[i]) {
        _X[i] = prevNonEmptyX;
        /* (((lastNonEmptyX >> 1) + lastNonEmptyC)
                 << 1); */  // we want to point to it's last element
      } else {
        prevNonEmptyX = _X[i];
      }
    }
    
  }

  // This select assumes an unoverse of _u/256 + m with no repeated elements
  // it returns what the _X[sq >> 8] would return
  // and the number of elements in the bucket
  pair<int64_t, uint32_t> select1(int64_t sq) const {
    // the position in _Y of the first byte of our relevant hint
    if (sq < 0) {
      return {-1, 0};
    }
    if (sq >= _n) {
      uint64_t last_el = _u;
      return {last_el, 0};
    }
    uint64_t h_idx = (sq / 32) * 36;
    uint64_t el_idx = h_idx + 4 + (sq % 32);
    uint64_t el = _Y[el_idx];
    uint64_t next_el_idx = ((sq + 1) / 32) * 36 + 4 + ((sq + 1) % 32);
    uint64_t next_el = _Y[next_el_idx];
    // uint32_t hint = ((uint32_t*)((_Y.data() + h_idx)))[0];
    uint32_t hint = *reinterpret_cast<const uint32_t*>(_Y.data() + h_idx);
    
    // hint is the upper bits of the element with rank sq
    // x is the number of elements before block _hint_
    uint32_t x = _X[hint];
    // non-empty bucket
    uint64_t y = x;  // this is where the rank of the first element of this
                     // bucket, that is how many ones before this
    uint64_t j = 1;
    uint64_t next_y = _X[hint + j];

    // we need to know when to stop scanning block so next_y has to point to
    // strictly after sq
    while (next_y <= sq && (hint + j) < _nblocks) {
      j++;
      next_y = _X[hint + j];
    }
    // now we scan from the upper bucket of our current element to find the
    // upper bucket of the next element
    uint64_t j_2 = j - 1;  // this is the bucket of the curr el
    uint64_t next_y_2 = _X[hint + j_2];
    while (next_y_2 <= (sq + 1) && (hint + j_2) < _nblocks) {
      j_2++;
      next_y_2 = _X[hint + j_2];
    }
    uint64_t upperbits = (hint + j - 1) << 8;
    uint64_t next_upperbits = ((hint + j_2 - 1) << 8);

    uint64_t curr_pos = upperbits + (uint64_t)el;
    uint64_t next_pos = next_upperbits + (uint64_t)next_el;
    // bcount will never get below zero as there are no repeated elements in our
    // set
    assert(next_pos > curr_pos);
    uint64_t bcount =
        (next_pos - curr_pos) -
        1;  // -1 because that is how many zeroes we had between them

    // update position to be returned as we don't want to include the count of
    // 0's.
    curr_pos -= sq;
    // if the block is empty, the position must point to the predecessor of the
    // block
    curr_pos -= (bool)(bcount == 0);

    return {(int64_t)curr_pos, (uint32_t)bcount};
  }

  size_t getu() const { return _u; }
  size_t getn() const { return _n; }

  uint64_t sizeInBytes() const {
    uint64_t sz = 3 * sizeof(uint64_t) + _X.size() * sizeof(uint32_t) +
                  _Y.size() * sizeof(uint8_t);
    return sz;
  }

  int64_t serialize(std::ostream& os) const {
    uint64_t written = 0;

    os.write(reinterpret_cast<const char*>(&_u), sizeof(_u));
    os.write(reinterpret_cast<const char*>(&_n), sizeof(_n));
    os.write(reinterpret_cast<const char*>(&_nblocks), sizeof(_nblocks));

    /* cout << "Pred8vS1_new::serialize: _X.size() = "
         << (_nblocks + 1) * sizeof(uint32_t) << " _Y.size() = " << _n << '\n'; */

    os.write(reinterpret_cast<const char*>(_X.data()),
             (_nblocks + 1) * sizeof(uint32_t));

    uint64_t n_hints = _n / 32 + 1;
    os.write(reinterpret_cast<const char*>(_Y.data()),
             (_n + n_hints * 4) * sizeof(uint8_t));

    written += 3 * sizeof(uint64_t) + sizeof(uint64_t) +
               _X.size() * sizeof(uint32_t) + _Y.size() * sizeof(uint8_t);

    return written;
  }

  void load(std::istream& is) {
    is.read(reinterpret_cast<char*>(&_u), sizeof(_u));
    is.read(reinterpret_cast<char*>(&_n), sizeof(_n));
    is.read(reinterpret_cast<char*>(&_nblocks), sizeof(_nblocks));

    _X.resize(_nblocks + 1);
    is.read(reinterpret_cast<char*>(_X.data()),
            (_nblocks + 1) * sizeof(uint32_t));
    uint64_t n_hints = _n / 32 + 1;
    _Y.resize(_n + n_hints * 4);  // +4 bytes per hint
    is.read(reinterpret_cast<char*>(_Y.data()),
            (_n + n_hints * 4) * sizeof(uint8_t));
  }

  Pred8vS1_new(Pred8vS1_new& other) {
    this->_u = other._u;
    this->_n = other._n;
    this->_nblocks = other._nblocks;
    this->_nActiveBuckets = other._nActiveBuckets;
    this->_X = other._X;
    this->_Y = other._Y;
  }

 private:
  uint64_t _u = 0;  // universe size
  uint64_t _n = 0;  // number of elements
  uint64_t _nblocks = 0;
  uint64_t _nActiveBuckets = 0;
  std::vector<uint32_t> _X;
  std::vector<uint8_t> _Y;
};

#endif
