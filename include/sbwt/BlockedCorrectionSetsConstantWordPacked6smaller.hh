#pragma once

#include <iostream>
#include <vector>

using namespace std;

template <int64_t sigma>
class BlockedCorrectionSetsBase4Rank3_ {
  uint64_t _logb = 6;
  uint64_t _b = (uint64_t)1 << _logb;  // number of symbols per block
  uint64_t _log_superb = 32;           // make this fit the blocksize
  uint64_t _super_b = (uint64_t)1 << _log_superb;
  uint64_t _n;
  uint64_t _N;
  uint64_t nsblocks;  // pointer just after the superblocks
  vector<uint64_t> _bits;
  vector<uint64_t> _bits_overflown;
  vector<uint32_t> _p;
  vector<uint8_t> _prefixsums;

  uint64_t n_blocks_in_superblock = _super_b / _b;
  uint64_t _log_ub = 16;
  uint64_t _ub = (uint64_t)1 << _log_ub;
  uint64_t n_blocks_in_ub = _ub / _b;
  uint64_t nublocks;

  uint64_t block_constant = 4 - 1;
  uint64_t CS_CONSTANT = 7;
  static constexpr size_t sums_per_prefix = 16;
  static constexpr size_t absoluteBytes = sigma * sizeof(uint16_t);
  static constexpr size_t groupBytes =
      absoluteBytes + sigma * (sums_per_prefix - 1);

 public:
  BlockedCorrectionSetsBase4Rank3_() {};

  BlockedCorrectionSetsBase4Rank3_(const std::string& seq,
                                   std::vector<uint64_t> correction_Set_A,
                                   std::vector<uint64_t> correction_Set_C,
                                   std::vector<uint64_t> correction_Set_G,
                                   std::vector<uint64_t> correction_Set_T,
                                   std::vector<uint64_t> correction_set_sizes) {
    if (sigma < 3 || sigma > 4) {
      std::cerr << "alphabet size = " << sigma << std::endl;
      throw std::invalid_argument("Works only for alphabets of size 3 or 4.");
    }
    _n = seq.size();
    uint64_t nblocks = (_n + _b - 1) / _b;
    nsblocks = _n / _super_b + 1;
    nublocks = _n / _ub + 1;
    uint64_t* block_sizes = new uint64_t[35];
    for (int k = 0; k < 35; k++) {
      block_sizes[k] = 0;
    }
    uint64_t* corr_totals = new uint64_t[_b + _b];
    for (int k = 0; k < _b + _b; k++) {
      corr_totals[k] = 0;
    }
    uint64_t fitted_block = 0;
    //_b is the number of subsets per block
    int64_t corrections = correction_Set_A.size() + correction_Set_C.size() +
                          correction_Set_G.size() + correction_Set_T.size();

    // length in bits of data structure; The +64 is for 4*16-bit ints per block
    // of _b 2-bit symbols. 64 bits for the lengths if the correction sets.
    // The last block has just the prefix sums (64 bits) but no other data.
    _N = nblocks * (block_constant * 64) + corrections * 16 + nblocks * 64;

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
    _p.reserve(nblocks + 2 + nsblocks * 8 + nblocks * 20);
    _p.resize(nblocks + 2 + nsblocks * 8 + nblocks * 20);
    size_t prefix_groups = (nblocks + sums_per_prefix - 1) / sums_per_prefix;

    _prefixsums.resize(prefix_groups * groupBytes);

    uint64_t* psums = new uint64_t[4];
    psums[0] = psums[1] = psums[2] = psums[3] = 0;
    uint64_t* super_sums = new uint64_t[4];
    super_sums[0] = super_sums[1] = super_sums[2] = super_sums[3] = 0;
    uint32_t* ub_sums = new uint32_t[4];
    ub_sums[0] = ub_sums[1] = ub_sums[2] = ub_sums[3] = 0;
    uint32_t* local_sums = new uint32_t[4];
    local_sums[0] = local_sums[1] = local_sums[2] = local_sums[3] = 0;

    uint64_t bi = 0;
    uint16_t* correction_set_lengths = new uint16_t[4];
    correction_set_lengths[0] = correction_set_lengths[1] =
        correction_set_lengths[2] = correction_set_lengths[3] = 0;

    uint64_t block_num = 0;
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
      // new block
      if (i % _b == 0) {
        /* cout << " Start of block, i = " << i << '\n'; */
        if (bi && bi - prev_size < min_words) {
          min_words = bi - prev_size;
        }
        if (bi - prev_size > max_words) {
          max_words = bi - prev_size;
        }
        block_sizes[bi - prev_size]++;
        bi = block_num * block_constant;

        prev_size = bi;

        // cout << "p_ptr " << p_ptr << "\n";

        for (int64_t cs = 0; cs < 4; cs++) {
          int64_t correction = (cs == 0 ? -correction_set_lengths[cs]
                                        : correction_set_lengths[cs]);
          psums[cs] += correction;
        }
        // Not storing prefix sums in bits data
        /* ((uint16_t*)(_bits.data() + bi))[0] = (uint16_t)psums[0];
        ((uint16_t*)(_bits.data() + bi))[1] = (uint16_t)psums[1];
        ((uint16_t*)(_bits.data() + bi))[2] = (uint16_t)psums[2];
        ((uint16_t*)(_bits.data() + bi))[3] = (uint16_t)psums[3];

        bi++;  // move past the words containing the
        // construct the correction set lengths for this block */

        // storing prefix sums in dedicated table for now
        // storing every 16th pref sum as 16 bits followed by
        // 15 8 bit prev block sums
        size_t pref_ind = block_num / sums_per_prefix;
        size_t inner_pref_ind = block_num % sums_per_prefix;
        if (block_num % sums_per_prefix == 0) {
          size_t groupBase = pref_ind * groupBytes;

          auto* pref16 =
              reinterpret_cast<uint16_t*>(_prefixsums.data() + groupBase);

          for (size_t s = 0; s < sigma; ++s) {
            pref16[s] = static_cast<uint16_t>(psums[s]);
            local_sums[s] = psums[s];
          }
        } else {
          size_t deltaBase = pref_ind * groupBytes + absoluteBytes;

          for (size_t s = 0; s < sigma; ++s) {
            auto delta = psums[s] - local_sums[s];
            assert(delta <= UINT8_MAX);

            _prefixsums[deltaBase + s * (sums_per_prefix - 1) + inner_pref_ind -
                        1] = static_cast<uint8_t>(delta);

            local_sums[s] = psums[s];
          }
        }

        uint64_t total_corrections = 0;
        for (int64_t cs = 0; cs < 4; cs++) {
          uint64_t cs_size = correction_set_sizes[(i / _b + 1) * 4 + cs] -
                             (correction_set_sizes[(i / _b) * 4 + cs]);
          //((uint8_t*)(_bits.data() + bi))[cs] = (uint8_t)cs_size;
          correction_set_lengths[cs] = cs_size;
          total_corrections += cs_size;
          if (cs_size == _b) {
            cout << "Correction set " << cs << " size: " << cs_size << " pos "
                 << i << '\n';
            cout << " next size " << correction_set_sizes[(i / _b + 1) * 4 + cs]
                 << " current size "
                 << (correction_set_sizes[(i / _b) * 4 + cs]) << "\n";
          }

          // Update psums
        }

        corr_totals[total_corrections]++;

        if (total_corrections > CS_CONSTANT) {
          // the first 32 bits: highest bit is 1. Then correction set sizes for
          // C, A, G, T, because C is the smallest
          ((uint8_t*)(_bits.data() + bi))[0] =
              (uint8_t)correction_set_lengths[1] | 1 << 7;
          ((uint8_t*)(_bits.data() + bi))[1] =
              (uint8_t)correction_set_lengths[0];
          ((uint8_t*)(_bits.data() + bi))[2] =
              (uint8_t)correction_set_lengths[2];
          ((uint8_t*)(_bits.data() + bi))[3] =
              (uint8_t)correction_set_lengths[3];

          // the next 32 bits is the pointer to the position in P array
          ((uint32_t*)(_bits.data() + bi))[1] = p_ptr;
          int64_t local_i = 0;
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
              ((uint8_t*)(_p.data() + p_ptr))[local_i] = bitpos;
              local_i++;
            }
          }
          bi++;
          p_ptr += (local_i + 3) / 4;

        } else {
          // the first 4 bits of the first byte is the number of correction
          // sets: at most 0111, so the highest bit = 0
          ((uint8_t*)(_bits.data() + bi))[0] = total_corrections << 4;
          // now we mark which correction sets are non-empty
          uint8_t cs_flags = 0x0;
          for (int64_t cs = 0; cs < 4; cs++) {
            if (correction_set_lengths[cs]) cs_flags |= (1 << cs);
          }
          // 1st byte: 0111 TGCA
          ((uint8_t*)(_bits.data() + bi))[0] |= cs_flags;

          int64_t local_i = 1;

          // Now the correction set positions ranges [0.._b) using _logb bits
          // We set the highest bit to indicate start of a new correction set
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
              if (!_i) {
                bitpos |= 1 << 7;
              }
              ((uint8_t*)(_bits.data() + bi))[local_i] = bitpos;
              local_i++;
            }
          }
          bi += (local_i + 7) / 8;  // move past the words containing the
          fitted_block++;
        }
        // correction sets for this block
        block_num++;
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
      //_p[p_ptr++] = (uint32_t)bi;
      bi++;
      cout << "  --- !!!!! this happened\n";
    }
    _bits.resize(bi + 1);
    _N = _bits.size() * 64;
    _p.resize(p_ptr + 1);
    std::cout << "Finished constructing BlockedCorrectionSetsBase4Rank3_"
              << " of size " << size_in_bytes() << " bytes" << std::endl;
    std::cout << "P array size: " << _p.size() * sizeof(uint32_t) << " bytes"
              << std::endl;
    std::cout << "_bits array size: " << _bits.size() * sizeof(uint64_t)
              << " bytes" << std::endl;
    std::cout << "_prefixsums array size: "
              << _prefixsums.size() * sizeof(uint8_t) << " bytes" << std::endl;

    cout << "Fitted " << fitted_block << " blocks, "
         << (float)fitted_block / block_num * 100 << "'%'  out of " << block_num
         << "\n";

    delete[] psums;
    delete[] correction_set_lengths;
    delete[] ub_sums;
    delete[] corr_totals;
    delete[] block_sizes;
    delete[] super_sums;
    delete[] local_sums;
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
    uint64_t block_num = pos >> _logb;

    uint64_t blockstart = block_num * block_constant;

    uint64_t super_sum =
        ((uint64_t*)(_p.data() + 8 * (pos >> _log_superb)))[sym];
    uint64_t ub_sum =
        ((uint32_t*)(_p.data() + nsblocks * 8 + 4 * ((pos >> _log_ub))))[sym];

    // Prefix-sum table layout (68 bytes per group of 16 blocks):
    //   bytes  0-7  : four uint16_t absolute prefix sums
    //   bytes  8-22 : 15 uint8_t deltas for symbol 0
    //   bytes 23-37 : 15 uint8_t deltas for symbol 1
    //   bytes 38-52 : 15 uint8_t deltas for symbol 2
    //   bytes 53-67 : 15 uint8_t deltas for symbol 3

    size_t group = block_num / sums_per_prefix;
    size_t offset = block_num % sums_per_prefix;

    const auto* pref16 = reinterpret_cast<const uint16_t*>(_prefixsums.data());

    size_t absIndex = group * (groupBytes / sizeof(uint16_t));

    // Absolute prefix sum.
    uint64_t preBlockRank = pref16[absIndex + sym];

    // Delta bytes.
    size_t deltaBase =
        group * groupBytes + absoluteBytes + sym * (sums_per_prefix - 1);

    for (size_t i = 0; i < offset; ++i) {
      preBlockRank += _prefixsums[deltaBase + i];
    }
    /* std::cout << "Querying rank for pos: " << pos << " symbol: " << sym
    << " preBlockRank: " << preBlockRank << '\n'; */

    int32_t correction = 0;

    bool jumped = 0;
    uint32_t corr_lens_a = 0;
    uint32_t corr_lens_c = 0;
    uint32_t corr_lens_g = 0;
    uint32_t corr_lens_t = 0;
    uint64_t correction_offset = 0;
    uint64_t correction_end = 0;

    uint8_t first_byte = ((uint8_t*)(_bits.data() + blockstart))[0];

    if (first_byte & 1 << 7) {
      jumped = 1;
      corr_lens_a = ((uint8_t*)(_bits.data() + blockstart))[1];
      corr_lens_c = first_byte - (1 << 7);
      corr_lens_g = ((uint8_t*)(_bits.data() + blockstart))[2];
      corr_lens_t = ((uint8_t*)(_bits.data() + blockstart))[3];
      correction_offset = ((sym > 0) * corr_lens_a + (sym > 1) * corr_lens_c +
                           (sym > 2) * corr_lens_g);
      correction_end = (corr_lens_a + (sym > 0) * corr_lens_c +
                        (sym > 1) * corr_lens_g + (sym > 2) * corr_lens_t);
    }

    uint32_t p_ptr;
    if (!jumped) {
      // retrieve the number total size of correction sets
      uint8_t total_corrections = first_byte >> 4;
      // get the lower 4 bits indicating which sets are present
      uint8_t corr_present = first_byte & 0xF;
      // cout << " First byte ";
      // print64bitword(first_byte);
      // get the correction
      if (corr_present & (1 << sym)) {
        // we actually have to do something
        uint8_t ones_before_me = 1;
        for (int8_t symi = 0; symi < sym; symi++) {
          ones_before_me += bool(corr_present & 1 << symi);
        }
        // cout << " I am present " << sym << " need to see ones: " <<
        // ones_before_me << "\n";
        uint8_t ones_seen = 0;
        for (uint64_t i = 0; i < total_corrections; i++) {
          uint16_t bitpos = ((uint8_t*)(_bits.data() + blockstart))[1 + i];
          if (bitpos & 1 << 7) ones_seen++;
          if (ones_before_me == ones_seen) {
            // get the 7 lowest bits
            bitpos &= 0x7F;
            if (bitpos >= (pos & (_b - 1))) {
              break;
            }
            correction++;
          }
          /* if (bitpos < (pos & (_b - 1))) {
            correction++;
            } else {
              break;
              } */
        }
      }

    } else {
      p_ptr = ((uint32_t*)(_bits.data() + blockstart))[1];
      // cout << "Looking for corrections elsewhere " << p_ptr << "\n";
      // cout << " correction_offset, correction_end, cs_size "
      //     << correction_offset << " " << correction_end << "\n";
      // get the correction
      //  __builtin_prefetch(&((uint8_t*)(_p.data() +
      //  p_ptr))[correction_offset]);
      /* for (uint64_t i = correction_offset; i < correction_end; i++) {
        uint16_t bitpos = ((uint8_t*)(_p.data() + p_ptr))[i];
        if (bitpos < (pos & (_b - 1))) {
          correction++;
        } else {
          break;
        }
      } */
    }

    /* const uint64_t* blockwords =
        _bits.data() + blockstart +
        (corr_lens_a + corr_lens_c + corr_lens_g + corr_lens_t + 4 + 7) /
            8;  // move past the prefix sums and correction lengths */
    // uint64_t added = !jumped * (total_corrections + 4 + 7) / 8 + jumped;
    /* uint32_t added = (total_corrections + 4 + 7) / 8;
    if (jumped) added = 1; */
    uint32_t added = 1;

    const uint64_t* blockwords =
        _bits.data() + blockstart +
        added;  // move past the prefix sums and correction lengths
    uint64_t blocki =
        ((pos & (_b - 1)) / 64) *
        2;  // index of word in this block containing the query position

    uint64_t leftOverRank = 0;

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

    if (jumped) {
      for (uint64_t i = correction_offset; i < correction_end; i++) {
        uint16_t bitpos = ((uint8_t*)(_p.data() + p_ptr))[i];
        if (bitpos < (pos & (_b - 1))) {
          correction++;
        } else {
          break;
        }
      }
    }

    int64_t result = preBlockRank + leftOverRank +
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
    uint64_t prefixsums_size = _prefixsums.size();
    out.write(reinterpret_cast<const char*>(&prefixsums_size),
              sizeof(uint64_t));

    out.write(reinterpret_cast<const char*>(_prefixsums.data()),
              sizeof(uint8_t) * _prefixsums.size());
    written += sizeof(uint64_t);
    written += sizeof(uint64_t);
    written += sizeof(uint32_t);
    written += sizeof(uint64_t) * (bits_size);
    written += sizeof(uint32_t) * (_p.size());
    written += sizeof(uint8_t) * _prefixsums.size();
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

    _p.resize(p_size);
    in.read((char*)_p.data(), sizeof(uint32_t) * (_p.size()));
    uint64_t prefixsums_size;
    in.read(reinterpret_cast<char*>(&prefixsums_size), sizeof(uint64_t));
    _prefixsums.resize(prefixsums_size);

    in.read(reinterpret_cast<char*>(_prefixsums.data()),
            sizeof(uint8_t) * _prefixsums.size());
    nsblocks = _n / _super_b + 1;
    nublocks = _n / _ub + 1;
    cout << "Constant size correction sets split, max corrections per block: "
         << CS_CONSTANT << " block size " << block_constant
         << " words, _b: " << _b << " \n";
  }
};
