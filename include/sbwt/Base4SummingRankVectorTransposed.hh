#pragma once

#include <iostream>
#include <vector>

using namespace std;

// For a sequence of lenghts from set {0,2,3,4}
// with mapping 0:0, 2:1, 3:2, 4:3
template <int64_t sigma>
class Base4SummingRankVectorTransposed {
  uint64_t _logb = 7;
  uint64_t _b = 1<<_logb;  // number of symbols per block
  uint64_t _n;
  uint64_t _N;
  vector<uint64_t> _bits;

 public:
  Base4SummingRankVectorTransposed() {};
  // Base4SummingRankVectorTransposed(const vector<char>& seq) {
  Base4SummingRankVectorTransposed(const std::string& seq) {
    if (sigma < 3 || sigma > 4) {
      std::cerr << "alphabet size = " << sigma << std::endl;
      throw std::invalid_argument("Works only for alphabets of size 3 or 4.");
    }
    _n = seq.size();
    uint64_t nblocks = (_n + _b - 1) / _b;
    //_b is the number of symbols per block

    // length in bits of data structure; The +64 is for 2*32-bit ints per block
    // of _b 2-bit symbols. The last block is only 1 word containing the pref sums.
    _N = nblocks * (2 * _b + 64) + 64;  // +64 for the last block's prefix sums

    cout << "_n: " << _n << " nblocks: " << nblocks << " _N: " << _N << '\n';
    _bits.reserve(_N / 64);
    _bits.resize(_N / 64);
    uint32_t psums = 0;
    uint64_t bi = 0;
    uint64_t psum_bi = bi;
    for (uint64_t i = 0; i < _n;) {
      // cout << "i = " << i << '\n';
      // keeping track of half block pref sums
      if (i % _b == 0) {
        ((uint32_t*)(_bits.data() + bi))[0] = psums;
        ((uint32_t*)(_bits.data() + bi))[1] = psums; // init now, update later
        // cout << "Setting prefix sum for 0 to: " << psums[0] << '\n';
        psum_bi = bi;
        bi++;  // move past the word containing the prefix sums for this block
      }
      if (i % _b && i % (_b >> 1) == 0) {
        ((uint32_t*)(_bits.data() + psum_bi))[1] = psums;
        // cout << "Setting prefix sum for 0 to: " << psums[0] << '\n';
      }

      uint64_t j = 0;
      uint64_t upper_w = 0;
      uint64_t lower_w = 0;
      int smalls = 0;
      int highs = 0;
      while (j < 64 && (i + j) < _n) {
        uint8_t sym = seq[i + j];
        psums += sym + (bool)sym;
        upper_w = upper_w | (((uint64_t)(bool)(sym & 0x2)) << (j));
        lower_w = lower_w | (((uint64_t)(bool)(sym & 0x1)) << (j));
        j++;
      }
      _bits[bi] = upper_w;
      bi++;
      _bits[bi] = lower_w;
      bi++;
      i += j;
    }

    if (_n % _b == 0) {
      // Set the psums of the last block
      uint32_t* last_block =
          (uint32_t*)(_bits.data() + nblocks * (2 * _b + 64) / 64);
      last_block[0] = psums;
    }
  }

  size_t size_in_bytes() const {
    size_t sz = 0;
    sz += (sizeof(uint64_t) * _bits.size());  //_bits
    sz += (sizeof(uint64_t) * 2);             //_b, _logb, _n, _N
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

  void summing_rank_in_payload(const uint64_t* blockwords, int64_t blocki,
                               uint64_t pos, uint64_t& wholeWordRank,
                               uint64_t& leftOverRank,
                               uint64_t start_pos = 0) const {
    for (uint64_t i = start_pos; i < blocki; i += 2) {
      uint64_t upper_w = blockwords[i];
      uint64_t lower_w = blockwords[i + 1];
      wholeWordRank += __builtin_popcountll(upper_w) *
                       2;  // all 3s and 4s add a count of 2 to the sum
      wholeWordRank += __builtin_popcountll(lower_w) *
                       2;  // all 2s and 4s add a count of 2 to the sum
      // now add the remaining counts of 3s
      uint64_t threes = upper_w & lower_w;
      wholeWordRank += __builtin_popcountll(
          upper_w ^ threes);  // all 3s add a count of 1 to the sum
    }

    if (pos % 64) {  // possibly inspect part of the next word
      uint64_t upper_w = blockwords[blocki];
      // compute an appropriate shift
      uint32_t shift =
          64 - (pos % 64);  // pos%64 is never 0 inside this if statement
      uint64_t lower_w = blockwords[blocki + 1];
      upper_w = upper_w << shift;
      lower_w = lower_w << shift;
      leftOverRank += __builtin_popcountll(upper_w) *
                      2;  // all 3s and 4s add a count of 2 to the sum
      leftOverRank += __builtin_popcountll(lower_w) *
                      2;  // all 2s and 4s add a count of 2 to the sum
      // now add the remaining counts of 3s
      uint64_t threes = upper_w & lower_w;
      leftOverRank += __builtin_popcountll(upper_w ^ threes);
    }
  }

  // Sum of ranks (weighted) in half-open interval [0..pos)
  int64_t sum_of_ranks(int64_t pos) const {
    uint64_t blockstart = (pos >> _logb) * (2 * _b + 64) / 64;
    // blockstart is the word offset of the start of the block containing
    // position i
    // 2*_b is the number of bits from symbols in a block, because there are _b
    // items per block and 2 bits per symbol
    // 128 = 4*32 is the number of bits needed for the preblock ranks for each
    // symbol
    bool first_half = ((pos % _b) < _b >> 1);
    uint64_t preBlockRank =
        ((uint32_t*)(_bits.data() +
                     blockstart))[!first_half];  // retrieve the appropriate
                                                 // preblock rank
    /* std::cout << "Querying rank for pos: " << pos << " symbol: " << sym
              << " preBlockRank: " << preBlockRank << '\n'; */
    const uint64_t* blockwords =
        _bits.data() + blockstart + 1 + (first_half ? 0 : 1) * _b / 64;
    uint64_t blocki =
        ((pos & ((_b >> 1) - 1)) / 64) *
        2;  // index of word in this block containing the query position

    uint64_t wholeWordRank = 0, leftOverRank = 0;

    summing_rank_in_payload(blockwords, blocki, pos, wholeWordRank,
                            leftOverRank);

    int64_t result = preBlockRank + wholeWordRank + leftOverRank;
    /* std::cout << "Preblock rank : " << preBlockRank
              << " Whole word rank : " << wholeWordRank
              << " Leftover rank in partial block: " << leftOverRank
              << " sym = " << sym << " result = " << result << '\n'; */
    return result;
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
    cout << "_b " << _b << endl;
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
    cout << "Loaded summing rank struct with block size _b " << _b << endl;
  }
};
