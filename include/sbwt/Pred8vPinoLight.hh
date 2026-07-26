#ifndef _PRED8_VPINO_L_H_
#define _PRED8_VPINO_L_H_

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

class Pred8vPinoLight {
 public:
  Pred8vPinoLight() {}
  Pred8vPinoLight(const vector<uint64_t>& data) {
    _n = data.size();
    _min = data[0];
    _u = data[data.size() - 1] - _min;
    _nblocks = _u / 256 + ((_u % 256) > 0);

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
    cerr << "Pred8vPinoLight: _u _n _min _nblocks _nActiveBuckets: " << _u
         << ' ' << _n << ' ' << _min << ' ' << _nblocks << ' '
         << _nActiveBuckets << '\n';

    // _Y = new uint8_t[_n];
    _Y.resize(_n);

    uint64_t yi = 0;
    for (uint64_t i = 0; i < _n;) {
      uint64_t v = data[i] - _min;
      uint64_t bcount = _X[v >> 8];
      _X[v >> 8] = yi;
      // (yi << 1) | (bcount > 0);  // LSB==1 indicates bucket is non-empty
      if (bcount) {
        _C[v >> 8] = (uint8_t)(bcount - 1);
        for (uint64_t j = 0; j < bcount; j++) {
          v = data[i] - _min;
          _Y[yi++] = v & 255;
          i++;
        }
      }
    }
    assert(yi == _n);
    _X[_nblocks] = yi; // pointing past the last element in data
    
    uint32_t nextNonEmptyX = _X[_nblocks - 1];
    for (uint64_t i = _nblocks - 1; i > 0; i--) {
      if (!_X[i]) {
        _X[i] = nextNonEmptyX;  // we want to point to it's last element
      } else {
        nextNonEmptyX = _X[i];
      }
    }
    // _X[_nblocks] = (((_X[_nblocks - 1]) + _C[_nblocks - 1]));
    
    delete[] _C;

    std::vector<uint64_t> block_indices;
    // the escape index as the Pred8 structure always starts from an element
    uint64_t bi = 0;
    block_indices.push_back(bi);
    for (uint64_t i = 0; i < _nblocks; i++) {
      uint32_t x = _X[i];
      uint64_t y = x;
      uint64_t bcount = (_X[1 + i] - y);
      bi += 1 + bcount;
      block_indices.push_back(bi);
    }
    //block_indices.push_back(bi+=1);
    _U = Pred8vS1_new(block_indices);

    cout << "returned from Pred8vS1 constructor" << endl;
    uint64_t wrongs = 0;
    for (uint64_t bi = 0; bi < _nblocks; bi++) {
      uint32_t x = _X[bi];
      uint64_t y = x;
      uint64_t bcount = (_X[1 + bi] - y);
      y = (bcount) ? y : y-1;;
      auto _y = _U.select1(bi);

      // the returned element should point to the first position of the block
      // or, if block is empty to its predecessor
      if (y != _y.first) {
        wrongs++;
        cout << "Pred8vPinoLight: ERROR: bi = " << bi << " y = " << y
             << " _y.first = " << _y.first << " block empty = " << !(bcount)
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
    }
    wrongs = 0;
    for (uint64_t key = _min; key <= _min + _u; ++key) {
      auto p1 = getPred(key);
      auto p2 = getPred_Select(key);

      if (p1 != p2) {
        cout << "Mismatch for key " << key << '\n';
        cout << "getPred      : (" << p1.first << ", " << p1.second << ")\n";
        cout << "getPredSelect: (" << p2.first << ", " << p2.second << ")\n";
        wrongs++;
      }
      if (wrongs > 20) exit(1);
    }
    cerr << "Pred8vPinoLight: sizeInBytes(): " << sizeInBytes() << '\n';
  }

  //  p is the index of the predecessor in the set
  //  bool is 1 if the value at index p is equal to key, else 0
  pair<int64_t, bool> inline getPred(int64_t key) const {
    if (key < _min) {
      return {-1, false};
    }
    key = key - _min;
    uint32_t x = _X[key >> 8];
    uint64_t bcount = (_X[1 + (key >> 8)]) - x;
    if (bcount) {
      // non-empty bucket
      uint64_t y = x;  // this is where, in _Y, the bucket elements are located
      uint64_t k = key & 255;
      for (uint64_t j = 0; j < bcount; j++) {
        if (_Y[y + j] >= k) {
          return {y + j - (1 - (_Y[y + j] == k)), (_Y[y + j] == k)};
        }
      }
      // getting here means that we did not find anything in the bucket that
      // was
      // >= k this means the last element of the bucket is the predecessor of
      // k, and it is not equal to k
      return {y + bcount - 1, false};
    }
    // the bucket that key belongs to is empty
    return {(x - 1), false};
  }

  pair<int64_t, bool> getPred_Select(int64_t key) const {
    if (key < _min) {
      return {-1, false};
    }
    key = key - _min;
    int64_t bi = key >> 8;
    auto p = _U.select1(bi);
    uint64_t y = p.first;
    int64_t bcount = p.second;
    if (bcount) {
      uint64_t k = key & 255;
      for (uint64_t j = 0; j < bcount; j++) {
        if (_Y[y + j] >= k) {
          return {y + j - (1 - (_Y[y + j] == k)), (_Y[y + j] == k)};
        }
      }
      // getting here means that we did not find anything in the bucket that
      // was
      // >= k this means the last element of the bucket is the predecessor of
      // k, and it is not equal to k
      return {y + bcount - 1, false};
    }
    return {y, false};
  }

  int64_t rank(int64_t pos) const {
    pair<int64_t, bool> r = getPred_Select(pos);
    return (r.first + 1) - ((uint64_t)(r.second));
  }

  size_t getu() const { return _u; }
  size_t getn() const { return _n; }

  uint64_t sizeInBytes() const {
    return sizeof(_u) + sizeof(_n) + sizeof(_min) + sizeof(_nblocks) +
           _Y.size() + _U.sizeInBytes();
  }

  int64_t serialize(std::ostream& os) const {
    uint64_t written = 0;

    os.write(reinterpret_cast<const char*>(&_u), sizeof(_u));
    os.write(reinterpret_cast<const char*>(&_n), sizeof(_n));
    os.write(reinterpret_cast<const char*>(&_min), sizeof(_min));
    os.write(reinterpret_cast<const char*>(&_nblocks), sizeof(_nblocks));

    cout << "Pred8vPinoLight::serialize: _Y.size() = " << _n
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

  ~Pred8vPinoLight() = default;

  Pred8vPinoLight(Pred8vPinoLight& other) {
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
  std::vector<uint32_t> _X;
  std::vector<uint8_t> _Y;
  Pred8vS1_new _U;
};

#endif
