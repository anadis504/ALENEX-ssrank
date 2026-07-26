#pragma once

#include <iostream>
#include <vector>

using namespace std;

template <int64_t sigma>
class BlockedCorrectionSetsBase4RankWordPackedByte {
  uint64_t _logb = 8;
  uint64_t _b = (uint64_t)1 << _logb;  // number of symbols per block
  uint64_t _log_superb = 32;           // make this fit the blocksize
  uint64_t _super_b = (uint64_t)1 << _log_superb;
  uint64_t _n;
  uint64_t _N;
  uint64_t nsblocks;  // pointer just after the superblocks
  vector<uint64_t> _bits;
  vector<uint32_t> _p;
  uint64_t n_blocks_in_superblock = _super_b / _b;
  uint64_t _log_ub = 16;
  uint64_t _ub = (uint64_t)1 << _log_ub;
  uint64_t n_blocks_in_ub = _ub / _b;
  uint64_t nublocks;
  uint64_t p_min = 10;

 public:
  BlockedCorrectionSetsBase4RankWordPackedByte() {};

  BlockedCorrectionSetsBase4RankWordPackedByte(
      const std::string& seq, std::vector<uint64_t> correction_Set_A,
      std::vector<uint64_t> correction_Set_C,
      std::vector<uint64_t> correction_Set_G,
      std::vector<uint64_t> correction_Set_T,
      std::vector<uint64_t> correction_set_sizes) {
    if (sigma < 3 || sigma > 4) {
      std::cerr << "alphabet size = " << sigma << std::endl;
      throw std::invalid_argument("Works only for alphabets of size 3 or 4.");
    }
    _n = seq.size();
    uint64_t nblocks = _n / _b + 1;
    nsblocks = _n / _super_b + 1;
    nublocks = _n / _ub + 1;

    //_b is the number of subsets per block
    int64_t corrections = correction_Set_A.size() + correction_Set_C.size() +
                          correction_Set_G.size() + correction_Set_T.size();

    // length in bits of data structure; The +128 is for 4*32-bit ints per block
    // of _b 2-bit symbols. 64 bits for the lengths if the correction sets.
    // The last block has just the prefix sums (128 bits) but no other data.
    _N = nblocks * (2 * _b + 128 + 64) + corrections * 16 + nblocks * 64;

    cout << "Packed P array + packed prefix sums blocked correction sets split "
            "byte packed, block size "
         << _b << " \n";
    cout << "_n: " << _n << " nblocks: " << nblocks << " _N: " << _N
         << " corrections: " << corrections << " nsblocks: " << nsblocks
         << '\n';
    cout << "log_b: " << _logb << " _b: " << _b
         << " log_superb: " << _log_superb << " super_b: " << _super_b << '\n';
    _bits.reserve(_N / 64);
    _bits.resize(_N / 64);
    _p.reserve(nblocks + 2 + nsblocks * 8);
    _p.resize(nblocks + 2 + nsblocks * 8);

    uint64_t* psums = new uint64_t[4];
    psums[0] = psums[1] = psums[2] = psums[3] = 0;
    uint64_t* super_sums = new uint64_t[4];
    super_sums[0] = super_sums[1] = super_sums[2] = super_sums[3] = 0;
    uint32_t* ub_sums = new uint32_t[4];
    ub_sums[0] = ub_sums[1] = ub_sums[2] = ub_sums[3] = 0;

    uint64_t bi = 0;
    uint16_t* correction_set_lengths = new uint16_t[4];
    correction_set_lengths[0] = correction_set_lengths[1] =
        correction_set_lengths[2] = correction_set_lengths[3] = 0;

    uint64_t block_num = 0;
    uint64_t i = 0;
    uint64_t p_ptr = nsblocks * 8 + nublocks * 4;
    uint64_t p_global_value = 0;
    uint64_t min_words = 600;
    uint64_t max_words = 0;
    uint64_t prev_size = 0;
    for (uint64_t i = 0; i < _n;) {
      if (i % _super_b == 0) {
        super_sums[0] += psums[0] - correction_set_lengths[0] + ub_sums[0];
        super_sums[1] += psums[1] + correction_set_lengths[1] + ub_sums[1];
        super_sums[2] += psums[2] + correction_set_lengths[2] + ub_sums[2];
        super_sums[3] += psums[3] + correction_set_lengths[3] + ub_sums[3];
        ((uint64_t*)(_p.data() + (i >> _log_superb) * 8))[0] = super_sums[0];
        ((uint64_t*)(_p.data() + (i >> _log_superb) * 8))[1] = super_sums[1];
        ((uint64_t*)(_p.data() + (i >> _log_superb) * 8))[2] = super_sums[2];
        ((uint64_t*)(_p.data() + (i >> _log_superb) * 8))[3] = super_sums[3];
        correction_set_lengths[0] = correction_set_lengths[1] =
            correction_set_lengths[2] = correction_set_lengths[3] = 0;
        psums[0] = psums[1] = psums[2] = psums[3] = 0;
        ub_sums[0] = ub_sums[1] = ub_sums[2] = ub_sums[3] = 0;

        std::cout << "Processing block starting at position " << i
                  << " superblock_log " << _log_superb << '\n';
        cout << "Superblock prefix sums: " << super_sums[0] << " "
             << super_sums[1] << " " << super_sums[2] << " " << super_sums[3]
             << '\n';
      }
      if (i % _ub == 0) {
        ub_sums[0] += psums[0] - correction_set_lengths[0];
        ub_sums[1] += psums[1] + correction_set_lengths[1];
        ub_sums[2] += psums[2] + correction_set_lengths[2];
        ub_sums[3] += psums[3] + correction_set_lengths[3];
        ((uint32_t*)(_p.data() + nsblocks * 8 +
                     (block_num / n_blocks_in_ub) * 4))[0] = ub_sums[0];
        ((uint32_t*)(_p.data() + nsblocks * 8 +
                     (block_num / n_blocks_in_ub) * 4))[1] = ub_sums[1];
        ((uint32_t*)(_p.data() + nsblocks * 8 +
                     (block_num / n_blocks_in_ub) * 4))[2] = ub_sums[2];
        ((uint32_t*)(_p.data() + nsblocks * 8 +
                     (block_num / n_blocks_in_ub) * 4))[3] = ub_sums[3];

        psums[0] = psums[1] = psums[2] = psums[3] = 0;
        correction_set_lengths[0] = correction_set_lengths[1] =
            correction_set_lengths[2] = correction_set_lengths[3] = 0;
      }
      if (i % _b == 0) {
        /* cout << " Start of block, i = " << i << '\n'; */
        if (bi && bi - prev_size < min_words) {
          min_words = bi - prev_size;
        }
        if (bi - prev_size > max_words) {
          max_words = bi - prev_size;
        }
        prev_size = bi;
        if (block_num % 5 == 0) {
          _p[p_ptr] = (uint32_t)bi;
          p_global_value = (uint64_t)bi;
          p_ptr++;
        } else {
          uint64_t offset = (uint64_t)bi - p_global_value;
          uint64_t inner_i = block_num % 5;
          offset -= inner_i * p_min;

          if (offset > 0xFF) {
            cerr << "Error: offset " << offset << " too large for block_num "
                 << block_num << " inner_i: " << inner_i << '\n';
          }
          ((uint8_t*)(_p.data() + p_ptr))[inner_i - 1] = (uint8_t)(offset);
          if (inner_i == 4) {
            p_ptr += 1;
          }
        }
        // cout << "p_ptr " << p_ptr << "\n";

        for (int64_t cs = 0; cs < 4; cs++) {
          int64_t correction = (cs == 0 ? -correction_set_lengths[cs]
                                        : correction_set_lengths[cs]);
          psums[cs] += correction;
        }
        ((uint16_t*)(_bits.data() + bi))[0] = (uint16_t)psums[0];
        ((uint16_t*)(_bits.data() + bi))[1] = (uint16_t)psums[1];
        ((uint16_t*)(_bits.data() + bi))[2] = (uint16_t)psums[2];
        ((uint16_t*)(_bits.data() + bi))[3] = (uint16_t)psums[3];

        bi++;  // move past the words containing the
        // construct the correction set lengths for this block
        bool full_set = 0;
        for (int64_t cs = 0; cs < 4; cs++) {
          uint64_t cs_size = correction_set_sizes[(i / _b + 1) * 4 + cs] -
                             (correction_set_sizes[(i / _b) * 4 + cs]);
          ((uint8_t*)(_bits.data() + bi))[cs] = (uint8_t)cs_size;
          correction_set_lengths[cs] = cs_size;
          if (cs_size > 0xFF) {
            full_set = 1;
            cout << "Correction set " << cs << " size: " << cs_size << " pos "
                 << i << '\n';
          }

          // Update psums
        }
        if (full_set) {
          ((uint16_t*)(_bits.data() + bi))[0] = 0xFFFF;
          bi++;
          block_num++;
        } else {
          // bi++;  // move past the correction set lengths for this block
          //  construct the correction sets for this block
          int64_t local_i = 4;
          for (int64_t cs = 0; cs < 4; cs++) {
            uint64_t cs_size = correction_set_lengths[cs];
            uint64_t start_pos =
                correction_set_sizes[(i == 0 ? 0 : (i / _b) * 4) + cs];

            for (uint64_t _i = 0; _i < cs_size; _i++) {
              uint64_t pos = start_pos + _i;
              switch (cs) {
                case 0:
                  pos = correction_Set_A[start_pos + _i];
                  break;
                case 1:
                  pos = correction_Set_C[start_pos + _i];
                  break;
                case 2:
                  pos = correction_Set_G[start_pos + _i];
                  break;
                case 3:
                  pos = correction_Set_T[start_pos + _i];
                  break;
                default:
                  break;
              }
              uint8_t bitpos = pos % _b;
              ((uint8_t*)(_bits.data() + bi))[local_i] = bitpos;
              local_i++;
            }
          }
          bi += (local_i + 7) / 8;  // move past the words containing the
          // correction sets for this block
          block_num++;
        }
      }
      // Now construct the concat string bits
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
      for (int64_t cs = 0; cs < 4; cs++) {
        int64_t correction = (cs == 0 ? -correction_set_lengths[cs]
                                      : correction_set_lengths[cs]);
        psums[cs] += correction;
      }
      uint32_t* last_block = (uint32_t*)(_bits.data() + bi);
      last_block[0] = (uint16_t)psums[0];
      last_block[1] = (uint16_t)psums[1];
      last_block[2] = (uint16_t)psums[2];
      last_block[3] = (uint16_t)psums[3];
      // ((uint32_t*)(_p.data() + p_pointer))[(_n / _b) % 2] = bi << 1;
      _p[p_ptr++] = (uint32_t)bi;
      bi++;
    }
    _bits.resize(bi + 1);
    _N = _bits.size() * 64;
    _p.resize(p_ptr + 1);
    std::cout
        << "Finished constructing BlockedCorrectionSetsBase4RankWordPackedByte"
        << " of size " << size_in_bytes() << " bytes" << std::endl;
    std::cout << "P array size: " << _p.size() * sizeof(uint32_t) << " bytes"
              << std::endl;
    cout << "Min words per block: " << min_words
         << " Max words per block: " << max_words << '\n';
    delete[] psums;
    delete[] correction_set_lengths;
  }

  size_t size_in_bytes() const {
    size_t sz = 0;
    sz += (sizeof(uint64_t) * _bits.size());  //_bits
    sz += (sizeof(uint64_t) * 4);             //_b, _logb, _n, _N
    sz += (sizeof(uint64_t*));                //_bits's pointer
    sz += (sizeof(uint32_t) * _p.size());     //_p
    return sz;
  }

  // Useful during debugging
  void print64bitword(uint64_t w) const {
    for (uint64_t i = 0; i < 64; i++) {
      cout << ((w >> (63 - i)) & 1);
    }
    cout << '\n';
  }

  // Rank of symbol in half-open interval [0..pos)
  int64_t rank(int64_t pos, uint64_t sym) const {
    uint64_t super_sum =
        ((uint64_t*)(_p.data() + 8 * (pos >> _log_superb)))[sym];
    // uint64_t blockstart = _p[(pos >> _logb) + nsblocks * 8 + nublocks * 4];
    uint64_t ub_sum =
        ((uint32_t*)(_p.data() + nsblocks * 8 + 4 * ((pos >> _log_ub))))[sym];
    uint64_t block_num = pos >> _logb;
    uint64_t p_global_offset = block_num / 5 * 2 + nsblocks * 8 + nublocks * 4;
    uint64_t blockstart = _p[p_global_offset];
    if (block_num % 5) {
      uint64_t inner_i = block_num % 5;
      blockstart += ((uint8_t*)(_p.data() + p_global_offset + 1))[inner_i - 1] +
                    inner_i * p_min;
    }

    uint64_t preBlockRank = ((uint16_t*)(_bits.data() + blockstart))[sym];

    uint64_t correction = 0;
    uint64_t corr_lens_a = ((uint8_t*)(_bits.data() + blockstart + 1))[0];
    uint64_t corr_lens_c = ((uint8_t*)(_bits.data() + blockstart + 1))[1];
    uint64_t corr_lens_g = ((uint8_t*)(_bits.data() + blockstart + 1))[2];
    uint64_t corr_lens_t = ((uint8_t*)(_bits.data() + blockstart + 1))[3];

    if (((uint16_t*)(_bits.data() + blockstart + 1))[0] == 0xFFFF) {
      correction = sym ? 0 : (pos & (_b - 1));
      corr_lens_a = corr_lens_c = corr_lens_g = corr_lens_t = 0;

    } else {
      uint64_t correction_offset =
          ((sym > 0) * corr_lens_a + (sym > 1) * corr_lens_c +
           (sym > 2) * corr_lens_g);
      uint64_t correction_end =
          (corr_lens_a + (sym > 0) * corr_lens_c + (sym > 1) * corr_lens_g +
           (sym > 2) * corr_lens_t);

      // get the correction
      for (uint64_t i = correction_offset; i < correction_end; i++) {
        uint16_t bitpos = ((uint8_t*)(_bits.data() + blockstart + 1))[4 + i];
        if (bitpos < (pos & (_b - 1))) {
          correction++;
        } else {
          break;
        }
      }
    }

    const uint64_t* blockwords =
        _bits.data() + blockstart + 1 +
        (corr_lens_a + corr_lens_c + corr_lens_g + corr_lens_t + 4 + 7) /
            8;  // move past the prefix sums and correction lengths

    uint64_t blocki =
        ((pos & (_b - 1)) / 64) *
        2;  // index of word in this block containing the query position

    uint64_t wholeWordRank = 0, leftOverRank = 0;

    for (uint64_t i = 0; i < blocki; i += 2) {
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

    /* for (uint64_t i = 0; i < blocki; i += 2) {
      uint64_t upper_w = blockwords[i];
      uint64_t lower_w = blockwords[i + 1];
      uint64_t highs = __builtin_popcountll(upper_w);
      uint64_t lows = 64 - highs;

      if (sym < 2) {
        countpB = lows ? __builtin_popcountll(lower_w << highs) : 0;
        countpA = lows - countpB;

      } else {
        countpB = highs ? __builtin_popcountll(lower_w >> lows) : 0;
        countpA = highs - countpB;
      }
      wholeWordRank += (sym & 1) ? countpB : countpA;
    }

    if (pos % 64) {  // possibly inspect part of the next word
      uint64_t upper_w = blockwords[blocki];
      // compute an appropriate shift
      uint32_t shift =
          64 - (pos % 64);  // pos%64 is never 0 inside this if statement
      uint64_t lower_w = blockwords[blocki + 1];
      uint64_t highs = __builtin_popcountll(upper_w);
      uint64_t lows = 64 - highs;
      uint64_t needed_highs = __builtin_popcountll(upper_w << shift);
      uint64_t needed_lows = (pos % 64) - needed_highs;

      if (sym < 2) {
        countpB = needed_lows ? __builtin_popcountll((lower_w << highs)
                                                     << (lows - needed_lows))
                              : 0;  // shift to keep needed lows
        countpA = needed_lows - countpB;

      } else {
        countpB = needed_highs
                      ? __builtin_popcountll(((lower_w >> lows) << lows)
                                             << (highs - needed_highs))
                      : 0;  // shift to keep needed highs
        countpA = needed_highs - countpB;
      }
      leftOverRank = (sym & 1) ? countpB : countpA;
    } */
    int64_t result = preBlockRank + wholeWordRank + leftOverRank +
                     (sym == 0 ? -correction : correction) + super_sum + ub_sum;

    return result;
  }

  int64_t serialize(ostream& out) const {
    size_t written = 0;

    out.write((char*)&_n, sizeof(uint64_t));
    out.write((char*)&_N, sizeof(uint64_t));
    uint64_t bits_size = _bits.size();
    out.write((char*)&bits_size, sizeof(uint64_t));
    out.write((char*)_bits.data(), sizeof(uint64_t) * (bits_size));
    uint32_t p_size = _p.size();
    out.write((char*)&p_size, sizeof(uint32_t));
    out.write((char*)_p.data(), sizeof(uint32_t) * (_p.size()));
    written += sizeof(uint64_t);
    written += sizeof(uint64_t);
    written += sizeof(uint32_t);
    written += sizeof(uint64_t) * (bits_size);
    written += sizeof(uint32_t) * (_p.size());
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
    uint32_t p_size;
    in.read((char*)&p_size, sizeof(uint32_t));
    _p.reserve(p_size);
    _p.resize(p_size);
    in.read((char*)_p.data(), sizeof(uint32_t) * (_p.size()));
    nsblocks = _n / _super_b + 1;
    nublocks = _n / _ub + 1;
    cout << "Packed P array + packed prefix sums blocked correction sets split "
            "byte packed, block size "
         << _b << " \n";
  }
};
