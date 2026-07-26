#pragma once

#include <iostream>
#include <vector>

using namespace std;

template <int64_t sigma>
class Base4RankVector {
  uint64_t _logb = 10;
  uint64_t _b = 1 << _logb;  // number of symbols per block
  uint64_t _n;
  uint64_t _N;
  vector<uint64_t> _bits;
  uint64_t _MASK[4] = {0x5555555555555555, 0xAAAAAAAAAAAAAAAA, 0, 0};

 public:
  Base4RankVector() {};
  // Base4RankVector(const vector<char>& seq) {
  Base4RankVector(const std::string& seq) {
    if (sigma < 3 || sigma > 4) {
      std::cerr << "alphabet size = " << sigma << std::endl;
      throw std::invalid_argument("Works only for alphabets of size 3 or 4.");
    }
    _n = seq.size();
    uint64_t symsPerWord = 32;
    uint64_t nblocks = _n / _b + 1;
    //_b is the number of symbols per block

    // length in bits of data structure; The +128 is for 4*32-bit ints per block
    // of _b 2-bit symbols. The last block has just the prefix sums (128 bits)
    // but no other data.
    _N = nblocks * (2 * _b + 128);

    cout << "_n: " << _n << " nblocks: " << nblocks << " _N: " << _N << " bits.size(): " << _N / 64 << '\n';
    _bits.reserve(_N / 64);
    _bits.resize(_N / 64);
    uint32_t* psums = new uint32_t[4];
    psums[0] = psums[1] = psums[2] = psums[3] = 0;
    uint64_t bi = 0;
    for (uint64_t i = 0; i < _n;) {
      // cout << "i = " << i << '\n';
      if (i % _b == 0) {
        // cout << " Start of block, i = " << i << '\n';
        // we are at the start of a block - set its prefix sums
        ((uint32_t*)(_bits.data() + bi))[0] = psums[0];
        ((uint32_t*)(_bits.data() + bi))[1] = psums[1];
        ((uint32_t*)(_bits.data() + bi))[2] = psums[2];
        ((uint32_t*)(_bits.data() + bi))[3] = psums[3];
        // cout << "Setting prefix sum for 0 to: " << psums[0] << '\n';
        bi +=
            2;  // move past the words containing the prefix sums for this block
      }
      uint64_t j = 0;
      uint64_t w = 0;
      while (j < 32 && (i + j) < _n) {
        uint8_t sym = seq[i + j];
        if (sigma == 3)
          sym++;  // for ternary sequences, remap the alphabet from {0,1,2} to
                  // {1,2,3}
        psums[sym]++;
        w = w | (((uint64_t)sym) << (2 * j));
        j++;
      }
      _bits[bi] = w;
      bi++;
      i += j;
    }

    if (_n % _b == 0) {
      // Set the psums of the last block
      uint32_t* last_block =
          (uint32_t*)(_bits.data() + (nblocks - 1) * (2 * _b + 128) / 64);
      last_block[0] = psums[0];
      last_block[1] = psums[1];
      last_block[2] = psums[2];
      last_block[3] = psums[3];
    }
    delete[] psums;
  }

  size_t size_in_bytes() const {
    size_t sz = 0;
    sz += (sizeof(uint64_t) * _bits.size());  //_bits
    sz += (sizeof(uint64_t) * 4);             //_b, _logb, _n, _N
    sz += (sizeof(uint64_t*));                //_bits's pointer
    return sz;
  }

  // Useful during debugging
  void print64bitword(uint64_t w) const {
    for (uint64_t i = 0; i < 64; i++) {
      cout << ((w >> (63 - i)) & 1);
    }
    cout << '\n';
  }

  inline void count1110(uint64_t w, uint64_t& countpA,
                        uint64_t& countpB) const {
    uint64_t w2 = w & 0xAAAAAAAAAAAAAAAA;  // popcnt(w2) == # of 11s + # of 10s
    uint64_t w3 = (w2 >> 1);
    w3 = w3 & w;
    countpA = __builtin_popcountll(w3);
    countpB = __builtin_popcountll(w2) - countpA;
  }

  inline void count1110Partial(uint64_t w, uint64_t& countpA, uint64_t& countpB,
                               uint32_t shift) const {
    uint64_t w2 = w & 0xAAAAAAAAAAAAAAAA;  // popcnt(w2) == # of 11s + # of 10s
    uint64_t w3 = (w2 >> 1);
    w3 = w3 & w;
    countpA = __builtin_popcountll((w3 << shift) >> (shift - 2));
    countpB = __builtin_popcountll((w2 << shift) >> (shift - 2)) - countpA;
  }

  // Sum of ranks (weighted) in half-open interval [0..pos)
  int64_t sum_of_ranks(int64_t pos) const {
    uint64_t blockstart = (pos >> _logb) * (2 * _b + 128) / 64;
    // blockstart is the word offset of the start of the block containing
    // position i 2*_b is the number of bits from symbols in a block, because
    // there are _b items per block and 2 bits per symbol 128 = 4*32 is the
    // number of bits needed for the preblock ranks for each symbol
    uint64_t preBlockRankA = ((uint32_t*)(_bits.data() + blockstart))[0];
    uint64_t preBlockRankB = ((uint32_t*)(_bits.data() + blockstart))[1];
    uint64_t preBlockRankC = ((uint32_t*)(_bits.data() + blockstart))[2];
    uint64_t preBlockRankD = ((uint32_t*)(_bits.data() + blockstart))[3];

    const uint64_t* blockwords =
        _bits.data() + blockstart + 2;  // move past the prefix sums
    uint64_t blocki =
        (pos & (_b - 1)) /
        32;  // index of word in this block containing the query position

    uint64_t countpA = 0, countpB = 0, countpC = 0, countpD = 0;
    uint64_t countp[4] = {0, 0, 0, 0};
    for (uint64_t i = 0; i < blocki; i++) {
      uint64_t w = blockwords[i];
      count1110(w, countpD, countpC);
      w = ~w;
      count1110(w, countpA, countpB);
      countp[0] += countpA;
      countp[1] += countpB;
      countp[2] += countpC;
      countp[3] += countpD;
    }

    if (pos % 32) {  // possibly inspect part of the next word
      uint64_t w = blockwords[blocki];
      // compute an appropriate shift
      uint32_t shift =
          64 - 2 * (pos % 32);  // pos%32 is never 0 inside this if statement
      count1110Partial(w, countpD, countpC, shift);
      w = ~w;
      count1110Partial(w, countpA, countpB, shift);
      countp[0] += countpA;
      countp[1] += countpB;
      countp[2] += countpC;
      countp[3] += countpD;
    }
    // countp[0] += preBlockRankA;  // this is not needed for sum_of_ranks
    countp[1] += preBlockRankB;
    countp[2] += preBlockRankC;
    countp[3] += preBlockRankD;
    int64_t lens_sum = countp[1] * 2 + countp[2] * 3 + countp[3] * 4;
    return lens_sum;
  }

  // Rank of symbol in half-open interval [0..pos)
  int64_t rank(int64_t pos, uint64_t sym) const {
    uint64_t blockstart = (pos >> _logb) * (2 * _b + 128) / 64;
    // blockstart is the word offset of the start of the block containing
    // position i
    // 2*_b is the number of bits from symbols in a block, because there are _b
    // items per block and 2 bits per symbol
    // 128 = 4*32 is the number of bits needed for the preblock ranks for each
    // symbol
    uint64_t preBlockRank = ((
        uint32_t*)(_bits.data() +
                   blockstart))[sym];  // retrieve the appropriate preblock rank
    const uint64_t* blockwords =
        _bits.data() + blockstart + 2;  // move past the prefix sums
    uint64_t blocki =
        (pos & (_b - 1)) /
        32;  // index of word in this block containing the query position
    uint64_t countpA = 0, countpB = 0, wholeWordRank = 0, leftOverRank = 0;
    if (sym < 2) {
      for (uint64_t i = 0; i < blocki; i++) {
        uint64_t w = blockwords[i];
        w = ~w;
        count1110(w, countpA, countpB);
        uint64_t count = countpB;
        bool x = (sym == 0 || sym == 3);
        if (x) count = countpA;
        wholeWordRank += count;
      }
      if (pos % 32) {  // possibly inspect part of the next word
        uint64_t w = blockwords[blocki];
        // compute an appropriate shift
        uint32_t shift =
            64 - 2 * (pos % 32);  // pos%32 is never 0 inside this if statement
        w = ~w;
        count1110Partial(w, countpA, countpB, shift);
        uint64_t count = countpB;
        bool x = (sym == 0 || sym == 3);
        if (x) count = countpA;
        leftOverRank = count;
      }
    } else {  // sym >= 2
      for (uint64_t i = 0; i < blocki; i++) {
        uint64_t w = blockwords[i];
        count1110(w, countpA, countpB);
        uint64_t count = countpB;
        bool x = (sym == 0 || sym == 3);
        if (x) count = countpA;
        wholeWordRank += count;
      }
      if (pos % 32) {  // possibly inspect part of the next word
        uint64_t w = blockwords[blocki];
        // compute an appropriate shift
        uint32_t shift =
            64 - 2 * (pos % 32);  // pos%32 is never 0 inside this if statement
        count1110Partial(w, countpA, countpB, shift);
        uint64_t count = countpB;
        bool x = (sym == 0 || sym == 3);
        if (x) count = countpA;
        leftOverRank = count;
      }
    }
    return preBlockRank + wholeWordRank + leftOverRank;
  }

  int64_t serialize(ostream& out) const {
    size_t written = 0;

    out.write((char*)&_n, sizeof(uint64_t));
    out.write((char*)&_N, sizeof(uint64_t));
    uint64_t bits_size = _bits.size();
    out.write((char*)&bits_size, sizeof(uint64_t));
    out.write((char*)_bits.data(), sizeof(uint64_t) * (bits_size));
    written += sizeof(uint64_t);
    written += sizeof(uint64_t);
    written += sizeof(uint64_t);
    written += sizeof(uint64_t) * (bits_size);
    return written;
  }

  void load(istream& in) {
    in.read((char*)&_n, sizeof(uint64_t));
    in.read((char*)&_N, sizeof(uint64_t));
    uint64_t bits_size;
    in.read((char*)&bits_size, sizeof(uint64_t));
    _bits.reserve(bits_size);
    _bits.resize(bits_size);
    in.read((char*)_bits.data(), sizeof(uint64_t) * (bits_size));
  }
};
