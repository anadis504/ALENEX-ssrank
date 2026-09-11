#pragma once

#include <iostream>
#include <vector>

using namespace std;

template <int64_t sigma>
class Base4RankVectorTransposed {
  uint64_t _logb = 10;
  uint64_t _b = 1 << _logb;  // number of symbols per block
  uint64_t _n;
  uint64_t _N;
  vector<uint64_t> _bits;

 public:
  Base4RankVectorTransposed() {};
  // Base4RankVectorTransposed(const vector<char>& seq) {
  Base4RankVectorTransposed(const std::string& seq) {
    if (sigma < 3 || sigma > 4) {
      std::cerr << "alphabet size = " << sigma << std::endl;
      throw std::invalid_argument("Works only for alphabets of size 3 or 4.");
    }
    _n = seq.size();
    uint64_t nblocks = _n / _b + 1;
    //_b is the number of symbols per block

    // length in bits of data structure; The +128 is for 4*32-bit ints per block
    // of _b 2-bit symbols. The last block has just the prefix sums (128 bits)
    // but no other data.
    _N = nblocks * (2 * _b + 128);

    /* cout << "_n: " << _n << " nblocks: " << nblocks << " _N: " << _N << '\n'; */
    _bits.reserve(_N / 64);
    _bits.resize(_N / 64);
    uint32_t* psums = new uint32_t[4];
    psums[0] = psums[1] = psums[2] = psums[3] = 0;
    uint64_t bi = 0;
    for (uint64_t i = 0; i < _n;) {
      // cout << "i = " << i << '\n';
      if (i % _b == 0) {
        ((uint32_t*)(_bits.data() + bi))[0] = psums[0];
        ((uint32_t*)(_bits.data() + bi))[1] = psums[1];
        ((uint32_t*)(_bits.data() + bi))[2] = psums[2];
        ((uint32_t*)(_bits.data() + bi))[3] = psums[3];
        // cout << "Setting prefix sum for 0 to: " << psums[0] << '\n';
        bi +=
            2;  // move past the words containing the prefix sums for this block
      }

      uint64_t j = 0;
      uint64_t upper_w = 0;
      uint64_t lower_w = 0;
      int smalls = 0;
      int highs = 0;
      while (j < 64 && (i + j) < _n) {
        uint8_t sym = seq[i + j];
        psums[sym]++;
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

  void rank_in_payload(const uint64_t* blockwords, int64_t blocki, uint64_t pos,
                       uint64_t sym, uint64_t& wholeWordRank,
                       uint64_t& leftOverRank, uint64_t start_pos = 0) const {
    for (uint64_t i = start_pos; i < blocki; i += 2) {
      uint64_t upper_w = blockwords[i];
      uint64_t lower_w = blockwords[i + 1];
      upper_w = (sym & 0x2) ? upper_w : ~upper_w;
      lower_w = (sym & 0x1) ? lower_w : ~lower_w;
      wholeWordRank += __builtin_popcountll(upper_w & lower_w);
    }

    if (pos % 64) {  // possibly inspect part of the next word
      uint64_t upper_w = blockwords[blocki];
      // compute an appropriate shift
      uint32_t shift =
          64 - (pos % 64);  // pos%64 is never 0 inside this if statement
      uint64_t lower_w = blockwords[blocki + 1];
      upper_w = (sym & 0x2) ? upper_w : ~upper_w;
      lower_w = (sym & 0x1) ? lower_w : ~lower_w;
      leftOverRank = __builtin_popcountll((upper_w & lower_w) << shift);
    }
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
    /* std::cout << "Querying rank for pos: " << pos << " symbol: " << sym
              << " preBlockRank: " << preBlockRank << '\n'; */
    const uint64_t* blockwords =
        _bits.data() + blockstart + 2;  // move past the prefix sums
    uint64_t blocki =
        ((pos & (_b - 1)) / 64) *
        2;  // index of word in this block containing the query position

    uint64_t wholeWordRank = 0, leftOverRank = 0;

    rank_in_payload(blockwords, blocki, pos, sym, wholeWordRank, leftOverRank);

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
    // cout << "_b " << _b << endl;
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
    // cout << "Loaded quad rank struct with block size _b " << _b << endl;
  }
};
