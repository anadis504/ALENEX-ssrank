#pragma once

#include <iostream>
#include <vector>

using namespace std;

template <int64_t sigma>
class BlockedSplitBase4RankWordPackedWT {
  uint64_t _logb = 9;
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
  uint64_t p_min = 18;

 public:
  BlockedSplitBase4RankWordPackedWT() {};

  BlockedSplitBase4RankWordPackedWT(
      const std::string& seq, std::vector<uint64_t> non_singeltons_positions,
      std::vector<uint8_t> non_singleton_cols,
      std::vector<uint64_t> n_non_singeltons, uint64_t n, uint64_t seq_size) {
    if (sigma < 3 || sigma > 4) {
      std::cerr << "alphabet size = " << sigma << std::endl;
      throw std::invalid_argument("Works only for alphabets of size 3 or 4.");
    }
    _n = n;
    uint64_t nblocks = _n / _b + 1;
    nsblocks = _n / _super_b + 1;
    nublocks = _n / _ub + 1;

    //_b is the number of subsets per block
    int64_t n_b = n_non_singeltons[n_non_singeltons.size() - 1];

    // length in bits of data structure; The +128 is for 4*32-bit ints per block
    // of _b 2-bit symbols. 64 bits for the lengths if the correction sets.
    // The last block has just the prefix sums (128 bits) but no other data.
    _N = nblocks * (2 * _b + 128 + 64) + n_b * 16 + nblocks * 64;
    cerr << "_n: " << _n << " nblocks: " << nblocks << " _N: " << _N
         << " nonsingelton sets: " << n_b << " nsblocks: " << nsblocks << '\n';
    cerr << "log_b: " << _logb << " _b: " << _b
         << " log_superb: " << _log_superb << " super_b: " << _super_b << '\n';
    _bits.reserve(_N / 64);
    _bits.resize(_N / 64);
    _p.reserve(nblocks + 2 + nsblocks * 8 + nublocks * 4);
    _p.resize(nblocks + 2 + nsblocks * 8 + nublocks * 4);

    uint64_t* psums = new uint64_t[4];
    psums[0] = psums[1] = psums[2] = psums[3] = 0;
    uint64_t* super_sums = new uint64_t[4];
    super_sums[0] = super_sums[1] = super_sums[2] = super_sums[3] = 0;
    uint32_t* ub_sums = new uint32_t[4];
    ub_sums[0] = ub_sums[1] = ub_sums[2] = ub_sums[3] = 0;

    uint64_t bi = 0;
    uint64_t block_num = 0;
    uint64_t str_ptr = 0;
    uint64_t block_m = 0;  // number of nonsingelton in the current block
    uint64_t block_s = 0;  // number of singeltons in the current block
    uint64_t s_counter = 0;
    uint64_t i = 0;
    uint64_t p_ptr = nsblocks * 8 + nublocks * 4;
    uint64_t p_global_value = 0;
    uint64_t min_words = 600;
    uint64_t max_words = 0;
    uint64_t prev_size = 0;
    for (uint64_t sets_processed = 0; sets_processed < _n;) {
      if (sets_processed % _super_b == 0) {
        super_sums[0] += psums[0] + ub_sums[0];
        super_sums[1] += psums[1] + ub_sums[1];
        super_sums[2] += psums[2] + ub_sums[2];
        super_sums[3] += psums[3] + ub_sums[3];
        ((uint64_t*)(_p.data() + (block_num / n_blocks_in_superblock) * 8))[0] =
            super_sums[0];
        ((uint64_t*)(_p.data() + (block_num / n_blocks_in_superblock) * 8))[1] =
            super_sums[1];
        ((uint64_t*)(_p.data() + (block_num / n_blocks_in_superblock) * 8))[2] =
            super_sums[2];
        ((uint64_t*)(_p.data() + (block_num / n_blocks_in_superblock) * 8))[3] =
            super_sums[3];

        psums[0] = psums[1] = psums[2] = psums[3] = 0;
        ub_sums[0] = ub_sums[1] = ub_sums[2] = ub_sums[3] = 0;

        std::cout << "Processing block starting at position " << i
                  << " set_position " << sets_processed << '\n';
        cout << "Superblock prefix sums: " << super_sums[0] << " "
             << super_sums[1] << " " << super_sums[2] << " " << super_sums[3]
             << '\n';
      }
      if (sets_processed % _ub == 0) {
        ub_sums[0] += psums[0];
        ub_sums[1] += psums[1];
        ub_sums[2] += psums[2];
        ub_sums[3] += psums[3];
        ((uint32_t*)(_p.data() + nsblocks * 8 +
                     (block_num / n_blocks_in_ub) * 4))[0] = ub_sums[0];
        ((uint32_t*)(_p.data() + nsblocks * 8 +
                     (block_num / n_blocks_in_ub) * 4))[1] = ub_sums[1];
        ((uint32_t*)(_p.data() + nsblocks * 8 +
                     (block_num / n_blocks_in_ub) * 4))[2] = ub_sums[2];
        ((uint32_t*)(_p.data() + nsblocks * 8 +
                     (block_num / n_blocks_in_ub) * 4))[3] = ub_sums[3];

        psums[0] = psums[1] = psums[2] = psums[3] = 0;
      }
      if (sets_processed % _b == 0) {
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

        ((uint16_t*)(_bits.data() + bi))[0] = (uint16_t)psums[0];
        ((uint16_t*)(_bits.data() + bi))[1] = (uint16_t)psums[1];
        ((uint16_t*)(_bits.data() + bi))[2] = (uint16_t)psums[2];
        ((uint16_t*)(_bits.data() + bi))[3] = (uint16_t)psums[3];

        bi++;  // move past the words containing the

        int64_t local_i = 0;
        block_m = n_non_singeltons[block_num + 1] - n_non_singeltons[block_num];

        uint64_t start_pos = block_m ? n_non_singeltons[block_num] : 0;
        // write the 4 first bytes as follows:
        // 10 bits: block_m, next 10 (could be 9) bits: first position, last 4
        // bits: first column, then 8 bits for the size of this encoding: number
        // of words before the singeltons encoding starts

        uint64_t curr_pos = block_m ? non_singeltons_positions[start_pos] : 0;
        uint16_t prev_pos = curr_pos % _b;
        uint8_t curr_col = block_m ? non_singleton_cols[start_pos] : 0;

        for (int64_t cs = 0; cs < 4; cs++) {
          psums[cs] += (bool)(curr_col & (1 << cs));
        }

        ((uint8_t*)(_bits.data() + bi))[0] = (uint8_t)(block_m >> 2);
        ((uint8_t*)(_bits.data() + bi))[1] =
            (uint8_t)(block_m) << 6 | (prev_pos >> 4);
        ((uint8_t*)(_bits.data() + bi))[2] =
            (uint8_t)(prev_pos << 4) | curr_col;

        uint64_t bigs_count = 0;
        uint64_t giants_count = 0;
        uint64_t delta = 0;
        local_i += 4;

        for (uint64_t _i = start_pos + 1; _i < start_pos + block_m; _i++) {
          uint64_t pos = non_singeltons_positions[_i];

          uint16_t bitpos = pos % _b;
          delta = bitpos - prev_pos;
          prev_pos = bitpos;
          uint8_t col = non_singleton_cols[_i];

          if (delta > 0xF) {
            ((uint8_t*)(_bits.data() + bi))[local_i] = 0x0 | (uint8_t)col;
            local_i++;
            if (delta > 0xFF) {
              ((uint8_t*)(_bits.data() + bi))[local_i] = (uint8_t)(delta >> 8);
              local_i++;
              ((uint8_t*)(_bits.data() + bi))[local_i] =
                  (uint8_t)(delta & 0xFF);
              giants_count++;
            } else {
              ((uint8_t*)(_bits.data() + bi))[local_i] = (uint8_t)delta;
              bigs_count++;
            }
          } else {
            ((uint8_t*)(_bits.data() + bi))[local_i] =
                ((uint8_t)delta << 4) | (uint8_t)col;
          }
          for (int64_t cs = 0; cs < 4; cs++) {
            psums[cs] += (bool)(col & (1 << cs));
          }
          local_i++;
        }
        uint64_t encoding_size = (local_i + 7) / 8;
        uint64_t difference = encoding_size - (block_m / 8);
        ((uint8_t*)(_bits.data() + bi))[3] =
            (uint8_t)difference;  // store the size of the correction set

        bi += (local_i + 7) / 8;  // move past the words containing the
        // correction sets for this block

        block_s = _b - block_m;
        block_num++;
        s_counter = 0;
        sets_processed += block_m;
      }

      // Now construct the concat string bits
      uint64_t j = 0;
      uint64_t upper_w = 0;
      uint64_t lower_w = 0;
      while (j < 64 && (i + j) < seq_size && s_counter < block_s) {
        uint8_t sym = seq[str_ptr++];
        psums[sym]++;
        upper_w = upper_w | (((uint64_t)(bool)(sym & 0x2)) << (j));
        lower_w = lower_w | (((uint64_t)(bool)(sym & 0x1)) << (j));
        j++;
        s_counter++;
      }
      _bits[bi] = upper_w;
      bi++;
      _bits[bi] = lower_w;
      bi++;
      i += j;
      sets_processed += j;
    }

    if (_n % _b == 0) {
      // Set the psums of the last block
      uint32_t* last_block = (uint32_t*)(_bits.data() + bi);
      last_block[0] = (uint32_t)psums[0];
      last_block[1] = (uint32_t)psums[1];
      last_block[2] = (uint32_t)psums[2];
      last_block[3] = (uint32_t)psums[3];
      _p[p_ptr++] = (uint32_t)bi;
      bi += 2;
    }
    _bits.resize(bi + 1);
    _N = _bits.size() * 64;
    _p.resize(p_ptr + 1);  // trim the unused part of the _p array
    delete[] psums;
    delete[] super_sums;
    delete[] ub_sums;
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

  uint64_t get_block_start(int64_t pos) const {
    uint64_t block_num = pos >> _logb;
    uint64_t p_global_offset = block_num / 5 * 2 + nsblocks * 8 + nublocks * 4;
    uint64_t blockstart = _p[p_global_offset];
    if (block_num % 5) {
      uint64_t inner_i = block_num % 5;
      blockstart += ((uint8_t*)(_p.data() + p_global_offset + 1))[inner_i - 1] +
                    inner_i * p_min;
    }
    return blockstart;
  }

  void rank_in_payload(const uint64_t* blockwords, int64_t blocki, uint64_t pos,
                       uint8_t sym, uint64_t& wholeWordRank,
                       uint64_t& leftOverRank, uint64_t start_pos = 0) const {
    uint64_t countpA = 0, countpB = 0;
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
    uint64_t super_sum =
        ((uint64_t*)(_p.data() + 8 * (pos >> _log_superb)))[sym];
    // uint64_t blockstart = _p[(pos >> _logb) + nsblocks * 8 + nublocks * 4];
    uint64_t ub_sum =
        ((uint32_t*)(_p.data() + nsblocks * 8 + 4 * ((pos >> _log_ub))))[sym];
    uint64_t blockstart = get_block_start(pos);

    uint64_t preBlockRank = ((
        uint16_t*)(_bits.data() +
                   blockstart))[sym];  // retrieve the appropriate preblock rank

    uint8_t* meta_ptr = ((uint8_t*)(_bits.data() + blockstart + 1));
    uint16_t block_m = meta_ptr[0];
    block_m = (uint16_t)(block_m << 2) | (meta_ptr[1] >> 6);
    uint16_t first_pos =
        (uint16_t)((meta_ptr[1] & 0x3F) << 4 | (meta_ptr[2] >> 4));
    uint8_t first_col = meta_ptr[2] & 0xF;
    uint8_t encoding_diff = meta_ptr[3] & 0x3F;
    uint32_t encoding_size = block_m / 8 + encoding_diff;

    // get the correction
    uint64_t sym_rank = 0;
    uint64_t non_singeltons_before_pos = 0;

    if (block_m && (first_pos < (pos & (_b - 1)))) {
      sym_rank += (bool)(first_col & (1 << sym));
      non_singeltons_before_pos++;
      for (uint64_t i = 0; non_singeltons_before_pos < block_m; i++) {
        uint8_t bits = ((uint8_t*)(_bits.data() + blockstart + 1))[i + 4];
        uint64_t delta = (bits >> 4);
        if (delta == 0) {
          i++;
          delta = ((uint8_t*)(_bits.data() + blockstart + 1))[i + 4];
          if (delta < 0x10) {
            delta =
                ((uint8_t*)(_bits.data() + blockstart + 1))[i + 5] + (1 << 8);
            i++;
          }
        }
        first_pos += delta;
        if ((first_pos) < (pos & (_b - 1))) {
          sym_rank += (bool)(bits & (1 << sym));
          non_singeltons_before_pos++;
        } else {
          break;
        }
      }
    }

    const uint64_t* blockwords =
        _bits.data() + blockstart + 1 +
        encoding_size;  // move past the prefix sums and bu length

    // our pos has to be updated to reflect the fact that we have removed the
    // non-singeltons
    pos -= non_singeltons_before_pos;

    uint64_t blocki =
        ((pos & (_b - 1)) >> 6)
        << 1;  // index of word in this block containing the query position
    uint64_t wholeWordRank = 0, leftOverRank = 0;
    rank_in_payload(blockwords, blocki, pos, sym, wholeWordRank, leftOverRank);
  
    int64_t result = preBlockRank + wholeWordRank + leftOverRank + sym_rank +
                     super_sum + ub_sum;

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
    cout << "Blocked split for large blocks, word packed, delta encoded, "
            "delta-P array, packed prefix sums, block size: "
         << _b << "\n";
  }
};
