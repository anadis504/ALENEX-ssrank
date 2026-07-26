#pragma once

#include <map>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include <vector>

#include "Base4RankVector.hh"
#include "Base4RankVectorWordPacked.hh"
#include "globals.hh"

namespace sbwt {

using namespace std;

template <typename quadrank_structure_t, typename bitvector_t,
          typename rank_support_t, typename pred8_t, typename pred16_t>
class SubsetCorrectionSetsRank {
  /* X_bitvector_type nonsingleton_sets;
  X_bitvector_rank_type nonsingleton_sets_rs; */

 public:
  quadrank_structure_t concat;
  pred8_t correction_Set_A_pred;
  pred16_t correction_Set_C_pred;
  pred16_t correction_Set_G_pred;
  pred16_t correction_Set_T_pred;

  // Count of character c in subsets up to pos, not including pos
  int64_t rank(int64_t pos, char c) const {
    /* std::cout << "Rank called for pos " << pos << " and char " << (int)c
              << '\n'; */

    uint8_t c_coded;
    int64_t correction = 0;
    switch (c) {
      case 'A':
        c_coded = 0;
        correction = -correction_Set_A_pred.rank(pos);
        break;
      case 'C':
        c_coded = 1;
        correction = correction_Set_C_pred.rank(pos);
        break;
      case 'G':
        c_coded = 2;
        correction = correction_Set_G_pred.rank(pos);
        break;
      case 'T':
        c_coded = 3;
        correction = correction_Set_T_pred.rank(pos);
        break;
      default:
        cerr << "Error: Rank called with non-ACGT character: " << c << endl;
        exit(1);
    }
    int64_t pre_result = concat.rank(pos, c_coded);
    /* std::cout << "Final result for rank(" << pos << "," << (int)c << ") is "
              << result << '\n'; */
    return pre_result + correction;
  }

  int64_t rank_by_charidx(int64_t pos, int64_t char_idx) const {
    int64_t correction = 0;
    int64_t pre_result = concat.rank(pos, char_idx);
    switch (char_idx) {
      case 0:
        correction = -correction_Set_A_pred.rank(pos);
        break;
      case 1:
        correction = correction_Set_C_pred.rank(pos);
        break;
      case 2:
        correction = correction_Set_G_pred.rank(pos);
        break;
      case 3:
        correction = correction_Set_T_pred.rank(pos);
        break;
      default:
        cerr << "Error: Rank called with non-ACGT character: " << char_idx
             << endl;
        exit(1);
    }
    return pre_result + correction;
  }

  bool contains(int64_t pos, char c) const {
    // TODO: faster
    int64_t r1 = this->rank(pos, c);
    int64_t r2 = this->rank(pos + 1, c);
    return r1 != r2;
  }

  SubsetCorrectionSetsRank() {}

  SubsetCorrectionSetsRank(const sdsl::bit_vector& A_bits,
                           const sdsl::bit_vector& C_bits,
                           const sdsl::bit_vector& G_bits,
                           const sdsl::bit_vector& T_bits) {
    assert(A_bits.size() == C_bits.size() && C_bits.size() == G_bits.size() &&
           G_bits.size() == T_bits.size());

    int64_t n = A_bits.size();
    int64_t n_b = 0;  // Number of branching nodes plus the nodes that do not
                      // have outedges
    int64_t n_u = 0;  // Number of nodes with exactly one outgoing edge
    for (int64_t i = 0; i < n; i++) {
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 1)
        n_u++;
      else {
        n_b++;
      }
    }

    std::string Y_str(n, '\0');
    std::vector<uint64_t> correction_Set_A;
    std::vector<uint64_t> correction_Set_C;
    std::vector<uint64_t> correction_Set_G;
    std::vector<uint64_t> correction_Set_T;

    int64_t Y_str_idx = 0;
    for (int64_t i = 0; i < n; i++) {
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

    concat = quadrank_structure_t(Y_str);

    correction_Set_A_pred = pred8_t(correction_Set_A);
    correction_Set_C_pred = pred16_t(correction_Set_C);
    correction_Set_G_pred = pred16_t(correction_Set_G);
    correction_Set_T_pred = pred16_t(correction_Set_T);

    // For debugging: verify that ranks match
    rank_support_t A_bits_rs;
    rank_support_t C_bits_rs;
    rank_support_t G_bits_rs;
    rank_support_t T_bits_rs;

    sdsl::util::init_support(A_bits_rs, &(A_bits));
    sdsl::util::init_support(C_bits_rs, &(C_bits));
    sdsl::util::init_support(G_bits_rs, &(G_bits));
    sdsl::util::init_support(T_bits_rs, &(T_bits));
    uint64_t wrongs = 0;
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
      if (wrongs > 20) exit(1);
    }
  }

  int64_t serialize(ostream& os) const {
    int64_t tmp, written = 0;
    written += concat.serialize(os);
    std::cout << "Serialized Quadrank_string " << written << " bytes\n";
    tmp = correction_Set_A_pred.serialize(os);
    std::cout << "Serialized correction_Set_A_pred " << tmp << " bytes\n";
    written += tmp;
    tmp = correction_Set_C_pred.serialize(os);
    std::cout << "Serialized correction_Set_C_pred " << tmp << " bytes\n";
    written += tmp;
    tmp = correction_Set_G_pred.serialize(os);
    std::cout << "Serialized correction_Set_G_pred " << tmp << " bytes\n";
    written += tmp;
    tmp = correction_Set_T_pred.serialize(os);
    std::cout << "Serialized correction_Set_T_pred " << tmp << " bytes\n";
    written += tmp;
    return written;
  }

  void load(istream& is) {
    concat.load(is);
    correction_Set_A_pred.load(is);
    correction_Set_C_pred.load(is);
    correction_Set_G_pred.load(is);
    correction_Set_T_pred.load(is);
  }

  SubsetCorrectionSetsRank(const SubsetCorrectionSetsRank& other) {
    assert(&other != this);  // What on earth are you trying to do?
    operator=(other);
  }

  SubsetCorrectionSetsRank& operator=(const SubsetCorrectionSetsRank& other) {
    if (&other != this) {
      this->concat = other.concat;
      this->correction_Set_A_pred = other.correction_Set_A_pred;
      this->correction_Set_C_pred = other.correction_Set_C_pred;
      this->correction_Set_G_pred = other.correction_Set_G_pred;
      this->correction_Set_T_pred = other.correction_Set_T_pred;
      return *this;
    } else {
      return *this;  // Assignment to self -> do nothing.
    }
  }
};

}  // namespace sbwt