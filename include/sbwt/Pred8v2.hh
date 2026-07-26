#ifndef _PRED8_V2_H_
#define _PRED8_V2_H_

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

class Pred8v2 {
   public:
    Pred8v2() {
    }
    Pred8v2(const vector<uint64_t> &data) {
       _n = data.size();
       _min = data[0];
       _u = data[data.size()-1] - _min;
       _nblocks = _u/256 + ((_u%256) > 1);

       _X = new uint32_t[_nblocks+1]; //+1 for a useful dummy at the end
       for(uint64_t i=0;i<_nblocks;i++) _X[i] = 0;
       uint8_t *_C = new uint8_t[_nblocks];
       for(uint64_t i=0;i<_nblocks;i++) _C[i] = 0;

       //NB: note that because we substract _min from everything, the 0th bucket is non-empty

       _nActiveBuckets = 0;
       for(uint64_t i=0;i<_n;i++){
          uint64_t v = data[i] - _min;
          if(!_X[v>>8]) _nActiveBuckets++;
          _X[v>>8]++;
       }

       _Y = new uint8_t[_n];

       cerr << "Pred8v2: _u _n _min _nblocks _nActiveBuckets: "<<_u<<' '<<_n<<' '<<_min<<' '<<_nblocks<<' '<<_nActiveBuckets<<'\n';
       cerr << "Pred8v2: sizeInBytes(): "<<sizeInBytes()<<'\n';
       uint64_t yi = 0;
       for(uint64_t i=0;i<_n;){
          uint64_t v = data[i] - _min;
          uint64_t bcount = _X[v >> 8];
          _X[v >> 8] = (yi << 1) | (bcount > 0); //LSB==1 indicates if bucket is non-empty
          if(bcount){
             _C[v >> 8] = (uint8_t)(bcount-1);
             for(uint64_t j = 0; j < bcount; j++){
                v = data[i] - _min;
                _Y[yi++] = v & 255;
                i++;
             }
          }//else{
          //   cerr << "Should never happen: "<<i<<'\n';
          //}
       }
       uint32_t lastNonEmptyX = _X[0] & 0xFFFFFFFE;
       uint32_t lastNonEmptyC = _C[0];
       for(uint64_t i=1;i<_nblocks;i++){
          if(!_X[i]){
             _X[i] = (((lastNonEmptyX>>1) + lastNonEmptyC)<<1); //we want to point to it's last element
          }else{
             lastNonEmptyX = _X[i] & 0xFFFFFFFE;
             lastNonEmptyC = _C[i];
          }
       }
       _X[_nblocks] = (((_X[_nblocks-1] >> 1) + _C[_nblocks-1])<<1);
       cerr << "_X[_nblocks]: " << (_X[_nblocks]>>1) << '\n';
    }


    //  p is the index of the predecessor in the set
    //  bool is 1 if the value at index p is equal to key, else 0
    pair<int64_t, bool> inline getPred(int64_t key) const{
        if (key < _min) {
            return {-1, false};
        }
        key = key - _min;
        uint32_t x = _X[key>>8];
        if(x&1){
           //non-empty bucket
           uint64_t y = x >> 1; //this is where, in _Y, the bucket elements are located, starting with the # of items in the bucket
           uint64_t bcount = ((_X[1+(key>>8)]>>1) + (1 - (_X[1+(key>>8)]&1))) - y;//_C[key>>8] + 1;
           //uint64_t bcount = _C[key>>8] + 1;
           //cerr << "bcount: "<<bcount<<'\n';
           //cerr << "_C[key>>8]+1: "<<(_C[key>>8] + 1)<<'\n';
           uint64_t k = key & 255;
           for(uint64_t j=0; j < bcount; j++){
              if(_Y[y+j] >= k){
                 return {y+j-(1-(_Y[y+j] == k)),(_Y[y+j] == k)};
              }
           }
           //getting here means that we did not find anything in the bucket that was >= k
           //this means the last element of the bucket is the predecessor of k, and it is not equal to k
           return {y+bcount-1,false};
        }
        //the bucket that key belongs to is empty
        return {(x>>1), false};
    }

    int64_t rank(int64_t pos) const{
       pair<int64_t, bool> r = getPred(pos);
       return (r.first+1)-((uint64_t)(r.second));
    } 

    size_t getu() const { return _u; }
    size_t getn() const { return _n; }

    uint64_t sizeInBytes() const{
        uint64_t sz = 4*sizeof(uint64_t) + (sizeof(uint32_t)*(_nblocks+1)) + _n;
        return sz;
    }

    int64_t serialize(std::ostream& os) const{
        int64_t written = 0;
        os.write((char *)&_u, sizeof(uint64_t));
        os.write((char *)&_n, sizeof(uint64_t));
        os.write((char *)&_min, sizeof(uint64_t));
        os.write((char *)&_nblocks, sizeof(uint64_t));
        os.write((char *)_X,sizeof(uint32_t)*(_nblocks+1));
        //os.write((char *)_C,sizeof(uint8_t)*_nblocks);
        os.write((char *)&_nActiveBuckets, sizeof(uint64_t));
        os.write((char *)_Y,sizeof(uint8_t)*_n);
        written += 5*sizeof(uint64_t) + (sizeof(uint32_t)*(_nblocks+1)) /* + _nblocks */ + _n;
        return written;
    }

    void load(std::istream& is){
       is.read((char *)&_u, sizeof(uint64_t));
       is.read((char *)&_n, sizeof(uint64_t));
       is.read((char *)&_min, sizeof(uint64_t));
       is.read((char *)&_nblocks, sizeof(uint64_t));
       _X = new uint32_t[_nblocks+1];
       is.read((char *)_X, sizeof(uint32_t)*(_nblocks+1));
       //_C = new uint8_t[_nblocks];
       //is.read((char *)_C, sizeof(uint8_t)*_nblocks);
       is.read((char *)&_nActiveBuckets, sizeof(uint64_t));
       _Y = new uint8_t[_n];
       is.read((char *)_Y, sizeof(uint8_t)*_n);
       //_X[_nblocks] = (((_X[_nblocks-1] >> 1) + _C[_nblocks-1])<<1);
    }
    
    Pred8v2(Pred8v2 &other){
       this->_u = other._u;
       this->_n = other._n;
       this->_min = other._min;
       this->_nblocks = other._nblocks;
       this->_nActiveBuckets = other._nActiveBuckets;
       this->_X = other._X;
       this->_Y = other._Y;
    }

   private:

    uint64_t _u = 0;  // universe size
    uint64_t _n = 0;  // number of elements
    uint64_t _min = 0;  // value of the smallest element
    uint64_t _nblocks = 0;
    uint64_t _nActiveBuckets = 0;
    uint32_t *_X;
    uint8_t *_Y;
    //uint8_t *_C;
};

#endif
