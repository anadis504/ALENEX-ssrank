#pragma once

#include <map>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include <vector>

// The blocked correction set with fixed block size
// #include "BlockedCorrectionSetsConstant.hh"
// The blocked correction set with block taking 2^8 subsets per block
#include "BlockedCorrectionSetsFullyPackedByte.hh"
// The Fixed block size block correction set with more corrections packed into on word only
// #include "BlockedCorrectionSetsConstantWordPacked.hh"
#include "Pred16.hh"
#include "Pred8v2.hh"
#include "globals.hh"

namespace sbwt {

using namespace std;

template <typename bitvector_t, typename rank_support_t>
class SubsetBlockedCorrectionSetsRank {
  /* X_bitvector_type nonsingleton_sets;
  X_bitvector_rank_type nonsingleton_sets_rs; */

 public:
  uint64_t _logb = 8;
  uint64_t _b = (uint64_t)1 << _logb;  // number of symbols per block
  BlockedCorrectionSetsBase4RankWordPackedByte<4> the_data_structure;

  // Count of character c in subsets up to pos, not including pos
  int64_t rank(int64_t pos, char c) const {
    /* std::cout << "Rank called for pos " << pos << " and char " << (int)c
              << '\n'; */

    uint8_t c_coded;
    switch (c) {
      case 'A':
        c_coded = 0;
        break;
      case 'C':
        c_coded = 1;
        break;
      case 'G':
        c_coded = 2;
        break;
      case 'T':
        c_coded = 3;
        break;
      default:
        cerr << "Error: Rank called with non-ACGT character: " << c << endl;
        exit(1);
    }
    int64_t result = the_data_structure.rank(pos, c_coded);
    /* std::cout << "Final result for rank(" << pos << "," << (int)c << ") is "
              << result << '\n'; */
    return result;
  }

  int64_t rank_by_charidx(int64_t pos, int64_t char_idx) const {
    return the_data_structure.rank(pos, char_idx);
  }

  bool contains(int64_t pos, char c) const {
    // TODO: faster
    int64_t r1 = this->rank(pos, c);
    int64_t r2 = this->rank(pos + 1, c);
    return r1 != r2;
  }

  SubsetBlockedCorrectionSetsRank() {}

  SubsetBlockedCorrectionSetsRank(const sdsl::bit_vector& A_bits,
                                  const sdsl::bit_vector& C_bits,
                                  const sdsl::bit_vector& G_bits,
                                  const sdsl::bit_vector& T_bits) {
    assert(A_bits.size() == C_bits.size() && C_bits.size() == G_bits.size() &&
           G_bits.size() == T_bits.size());

    int64_t n = A_bits.size();
    int64_t n_b = 0;  // Number of branching nodes plus the nodes that do not
                      // have outedges
    int64_t n_u = 0;  // Number of nodes with exactly one outgoing edge
    int64_t n_e = 0;  // Number of nodes with empty outgoing sets
    for (int64_t i = 0; i < n; i++) {
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 1)
        n_u++;
      else {
        n_b++;
        if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 0) n_e++;
      }
    }

    std::string Y_str(n, '\0');
    std::vector<uint64_t> correction_Set_A;
    std::vector<uint64_t> correction_Set_C;
    std::vector<uint64_t> correction_Set_G;
    std::vector<uint64_t> correction_Set_T;
    std::vector<uint64_t> correction_set_sizes((n / _b + 2) * 4, 0);

    int64_t Y_str_idx = 0;
    int64_t block_idx = 0;
    for (int64_t i = 0; i < n; i++) {
      if (i % _b == 0) {
        correction_set_sizes[block_idx * 4 + 0] = correction_Set_A.size();
        correction_set_sizes[block_idx * 4 + 1] = correction_Set_C.size();
        correction_set_sizes[block_idx * 4 + 2] = correction_Set_G.size();
        correction_set_sizes[block_idx * 4 + 3] = correction_Set_T.size();
        block_idx++;
      }
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 1) {
        // One outgoing label
        if (A_bits[i] == 1)
          Y_str[Y_str_idx++] = 0;
        else if (C_bits[i] == 1)
          Y_str[Y_str_idx++] = 1;
        else if (G_bits[i] == 1)
          Y_str[Y_str_idx++] = 2;
        else if (T_bits[i] == 1)
          Y_str[Y_str_idx++] = 3;
      } else {
        // 0 or > 1 outgoing labels
        if (!(A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i])) {
          // Empty set: assigning it to A and inserting to correction set A
          correction_Set_A.push_back(i);
          Y_str[Y_str_idx++] = 0;
          /* if (i >= 226823169 && i <= 226823169 + 20)  {
            std::cout << "Empty set at position " << i << "
          correction_set_sizes[block_idx * 4 + 0] = " <<
          correction_set_sizes[block_idx * 4 + 0] << '\n';
          } */
        } else if (A_bits[i] == 1) {
          // > 1 outgoing labels, where lex smallest is A
          Y_str[Y_str_idx++] = 0;
          // Adding index to correction sets of other present characters
          if (C_bits[i] == 1) {
            correction_Set_C.push_back(i);
          }
          if (G_bits[i] == 1) {
            correction_Set_G.push_back(i);
          }
          if (T_bits[i] == 1) {
            correction_Set_T.push_back(i);
          }
        } else if (C_bits[i] == 1) {
          // > 1 outgoing labels, where lex smallest is C
          Y_str[Y_str_idx++] = 1;
          if (G_bits[i] == 1) {
            correction_Set_G.push_back(i);
          }
          if (T_bits[i] == 1) {
            correction_Set_T.push_back(i);
          }
        } else if (G_bits[i] == 1) {
          // > 1 outgoing labels, where lex smallest is G and T has to be
          // present
          Y_str[Y_str_idx++] = 2;

          if (T_bits[i] == 1) {
            correction_Set_T.push_back(i);
          }
        }
      }
    }
    correction_set_sizes[block_idx * 4 + 0] = correction_Set_A.size();
    correction_set_sizes[block_idx * 4 + 1] = correction_Set_C.size();
    correction_set_sizes[block_idx * 4 + 2] = correction_Set_G.size();
    correction_set_sizes[block_idx * 4 + 3] = correction_Set_T.size();

    the_data_structure = BlockedCorrectionSetsBase4RankWordPackedByte<4>(
        Y_str, correction_Set_A, correction_Set_C, correction_Set_G,
        correction_Set_T, correction_set_sizes);

    // For debugging: verify that ranks match
    /*
    rank_support_t A_bits_rs;
    rank_support_t C_bits_rs;
    rank_support_t G_bits_rs;
    rank_support_t T_bits_rs;
    
    sdsl::util::init_support(A_bits_rs, &(A_bits));
    sdsl::util::init_support(C_bits_rs, &(C_bits));
    sdsl::util::init_support(G_bits_rs, &(G_bits));
    sdsl::util::init_support(T_bits_rs, &(T_bits));
    int wrongs = 0;
    for (int64_t i = 0; i < n; i++) {
      int64_t rA = A_bits_rs.rank(i);
      int64_t rC = C_bits_rs.rank(i);
      int64_t rG = G_bits_rs.rank(i);
      int64_t rT = T_bits_rs.rank(i);
      // std::cout << A_bits[i] << C_bits[i] << G_bits[i] << T_bits[i] << '\n';
      int64_t ownA = this->rank(i, 'A');
      int64_t ownC = this->rank(i, 'C');
      int64_t ownG = this->rank(i, 'G');
      int64_t ownT = this->rank(i, 'T');
      if (!(rA == ownA)) {
        std::cerr << "Rank mismatch at position " << i << " for A: " << rA
        << " vs " << ownA << '\n';
        wrongs++;
      }
      if (!(rC == ownC)) {
        std::cerr << "Rank mismatch at position " << i << " for C: " << rC
        << " vs " << ownC << '\n';
        wrongs++;
      }
      if (!(rG == ownG)) {
        std::cerr << "Rank mismatch at position " << i << " for G: " << rG
        << " vs " << ownG << '\n';
        wrongs++;
      }
      if (!(rT == ownT)) {
        std::cerr << "Rank mismatch at position " << i << " for T: " << rT
        << " vs " << ownT << '\n';
        wrongs++;
      }
      if (wrongs > 20) {
        std::cerr << "More than 20 wrong ranks, stopping checking.\n";
        exit(1);
        break;
      }
    } 
    */
  }

  int64_t serialize(ostream& os) const {
    int64_t tmp, written = 0;
    written += the_data_structure.serialize(os);
    std::cout << "Serialized the_data_structure " << written << " bytes\n";
    return written;
  }

  void load(istream& is) { the_data_structure.load(is); }

  SubsetBlockedCorrectionSetsRank(
      const SubsetBlockedCorrectionSetsRank& other) {
    assert(&other != this);  // What on earth are you trying to do?
    operator=(other);
  }

  SubsetBlockedCorrectionSetsRank& operator=(
      const SubsetBlockedCorrectionSetsRank& other) {
    if (&other != this) {
      this->the_data_structure = other.the_data_structure;
      return *this;
    } else {
      return *this;  // Assignment to self -> do nothing.
    }
  }
};

}  // namespace sbwt