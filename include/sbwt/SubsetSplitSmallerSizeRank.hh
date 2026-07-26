#pragma once

#include <map>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include <vector>

#include "Base4RankVector.hh"
#include "Base4RankVectorWordPacked.hh"
#include "Base4RankVectorTransposed.hh"
#include "Pred8v2.hh"
#include "globals.hh"

namespace sbwt {

using namespace std;

template <typename bitvector_t, typename rank_support_t>
class SubsetSplitSmallerSizeRank {
  /* X_bitvector_type nonsingleton_sets;
  X_bitvector_rank_type nonsingleton_sets_rs; */

 public:
  Base4RankVectorTransposed<4> concat;

  Pred8v2 nonsingleton_sets;
  /* Pred8v2 empty_sets; */

  bitvector_t empty_sets;
  rank_support_t empty_sets_rs;

  bitvector_t Z_A;
  bitvector_t Z_C;
  bitvector_t Z_G;
  bitvector_t Z_T;

  rank_support_t Z_A_rs;
  rank_support_t Z_C_rs;
  rank_support_t Z_G_rs;
  rank_support_t Z_T_rs;

  int64_t serialize(ostream& os) const {
    int64_t written = 0;
    int64_t written_tmp = 0;
    written += concat.serialize(os);
    cerr << "X serialized " << written << " bytes\n";
    written_tmp = nonsingleton_sets.serialize(os);
    cerr << " non-singelton sets vector serialized " << written_tmp << "\n";
    written += written_tmp;
    written_tmp = empty_sets.serialize(os);
    cerr << " empty sets vector serialized " << written_tmp << "\n";
    written += written_tmp;

    written_tmp = empty_sets_rs.serialize(os);
    cerr << " empty sets rs vector serialized " << written_tmp << "\n";
    written += written_tmp;

    written_tmp = Z_A.serialize(os);
    cerr << "Z_A serialized " << written_tmp << " bytes\n";
    written += written_tmp;
    written += Z_C.serialize(os);
    written += Z_G.serialize(os);
    written += Z_T.serialize(os);
    written += Z_A_rs.serialize(os);
    written += Z_C_rs.serialize(os);
    written += Z_G_rs.serialize(os);
    written += Z_T_rs.serialize(os);
    return written;
  }

  void load(istream& is) {
    concat.load(is);
    nonsingleton_sets.load(is);
    empty_sets.load(is);

    if (std::is_same<sdsl::rank_support_v5<>, rank_support_t>::value) {
      // Special handling needed for rank_support_v5 because of a design flaw in
      // sdsl
      empty_sets_rs.load(is, &empty_sets);
    } else {
      empty_sets_rs.load(is);
      empty_sets_rs.set_vector(&empty_sets);
    }
    Z_A.load(is);
    Z_C.load(is);
    Z_G.load(is);
    Z_T.load(is);

    if (std::is_same<sdsl::rank_support_v5<>, rank_support_t>::value) {
      // Special handling needed for rank_support_v5 because of a design flaw in
      // sdsl
      Z_A_rs.load(is, &Z_A);
      Z_C_rs.load(is, &Z_C);
      Z_G_rs.load(is, &Z_G);
      Z_T_rs.load(is, &Z_T);
    } else {
      Z_A_rs.load(is);
      Z_C_rs.load(is);
      Z_G_rs.load(is);
      Z_T_rs.load(is);
      Z_A_rs.set_vector(&Z_A);
      Z_C_rs.set_vector(&Z_C);
      Z_G_rs.set_vector(&Z_G);
      Z_T_rs.set_vector(&Z_T);
    }
  }

  SubsetSplitSmallerSizeRank(const SubsetSplitSmallerSizeRank& other) {
    assert(&other != this);  // What on earth are you trying to do?
    operator=(other);
  }

  SubsetSplitSmallerSizeRank& operator=(
      const SubsetSplitSmallerSizeRank& other) {
    if (&other != this) {
      this->concat = other.concat;
      this->nonsingleton_sets = other.nonsingleton_sets;
      this->empty_sets = other.empty_sets;
      this->empty_sets_rs = other.empty_sets_rs;
      this->empty_sets_rs.set_vector(&this->empty_sets);
      this->Z_A = other.Z_A;
      this->Z_C = other.Z_C;
      this->Z_G = other.Z_G;
      this->Z_T = other.Z_T;
      this->Z_A_rs = other.Z_A_rs;
      this->Z_C_rs = other.Z_C_rs;
      this->Z_G_rs = other.Z_G_rs;
      this->Z_T_rs = other.Z_T_rs;
      this->Z_A_rs.set_vector(&this->Z_A);
      this->Z_C_rs.set_vector(&this->Z_C);
      this->Z_G_rs.set_vector(&this->Z_G);
      this->Z_T_rs.set_vector(&this->Z_T);
      return *this;
    } else {
      return *this;  // Assignment to self -> do nothing.
    }
  }

  // Count of character c in subsets up to pos, not including pos
  int64_t rank(int64_t pos, char c) const {
    /* std::cout << "Rank called for pos " << pos << " and char " << (int)c
              << '\n'; */

    uint8_t c_coded;
    int64_t ns_rank = nonsingleton_sets.rank(pos);
    int64_t ne_rank = ns_rank - empty_sets_rs.rank(ns_rank);
    int64_t correction = 0;
    switch (c) {
      case 'A':
        c_coded = 0;
        correction = Z_A_rs.rank(ne_rank);
        break;
      case 'C':
        c_coded = 1;
        correction = Z_C_rs.rank(ne_rank);
        break;
      case 'G':
        c_coded = 2;
        correction = Z_G_rs.rank(ne_rank);
        break;
      case 'T':
        c_coded = 3;
        correction = Z_T_rs.rank(ne_rank);
        break;
      default:
        cerr << "Error: Rank called with non-ACGT character: " << c << endl;
        exit(1);
    }
    int64_t pre_result = concat.rank(pos - ns_rank, c_coded);
    /* std::cout << "Final result for rank(" << pos << "," << (int)c << ") is "
              << result << '\n'; */
    return pre_result + correction;
  }

  int64_t rank_by_charidx(int64_t pos, int64_t char_idx) const {
    int64_t ns_rank = nonsingleton_sets.rank(pos);
    int64_t ne_rank = ns_rank - empty_sets_rs.rank(ns_rank);
    int64_t correction = 0;
    switch (char_idx) {
      case 0:
        correction = Z_A_rs.rank(ne_rank);
        break;
      case 1:
        correction = Z_C_rs.rank(ne_rank);
        break;
      case 2:
        correction = Z_G_rs.rank(ne_rank);
        break;
      case 3:
        correction = Z_T_rs.rank(ne_rank);
        break;
      default:
        cerr << "Error: Rank called with non-ACGT character: " << char_idx
             << endl;
        exit(1);
    }
    int64_t pre_result =
        concat.rank(pos - ns_rank, char_idx_to_DNA(char_idx));
    /* std::cout << "Final result for rank(" << pos << "," << (int)c << ") is "
              << result << '\n'; */
    return pre_result + correction;
  }

  bool contains(int64_t pos, char c) const {
    // TODO: faster
    int64_t r1 = this->rank(pos, c);
    int64_t r2 = this->rank(pos + 1, c);
    return r1 != r2;
  }

  SubsetSplitSmallerSizeRank() {}

  SubsetSplitSmallerSizeRank(const sdsl::bit_vector& A_bits,
                             const sdsl::bit_vector& C_bits,
                             const sdsl::bit_vector& G_bits,
                             const sdsl::bit_vector& T_bits) {
    assert(A_bits.size() == C_bits.size() && C_bits.size() == G_bits.size() &&
           G_bits.size() == T_bits.size());

    int64_t n = A_bits.size();
    int64_t n_b = 0;  // Number of branching nodes plus the nodes that do not
                      // have outedges
    int64_t n_u = 0;  // Number of nodes with exactly one outgoing edge
    int64_t n_empty = 0;

    std::vector<uint64_t> nonsingleton_sets_vec;
    for (int64_t i = 0; i < n; i++) {
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 1)
        n_u++;
      else {
        nonsingleton_sets_vec.push_back(i);
        n_b++;
        if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 0) n_empty++;
      }
    }
    std::cout << "n: " << n << " n_b: " << n_b << " n_u: " << n_u
              << " n_empty: " << n_empty << '\n';

    std::string Y_str(n_u, '\0');
    sdsl::bit_vector empty_sets_plain(n_b);
    /* std::vector<uint64_t> empty_sets_plain; */
    sdsl::bit_vector Z_A_plain(n_b - n_empty);
    sdsl::bit_vector Z_C_plain(n_b - n_empty);
    sdsl::bit_vector Z_G_plain(n_b - n_empty);
    sdsl::bit_vector Z_T_plain(n_b - n_empty);

    int64_t Y_str_idx = 0;
    int64_t Z_idx = 0;
    int64_t nonsingletons_idx = 0;
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
        if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 0) {
          empty_sets_plain[nonsingletons_idx] = 1;
          /* empty_sets_plain.push_back(nonsingletons_idx); */
        } else {
          empty_sets_plain[nonsingletons_idx] = 0;
          Z_A_plain[Z_idx] = A_bits[i];
          Z_C_plain[Z_idx] = C_bits[i];
          Z_G_plain[Z_idx] = G_bits[i];
          Z_T_plain[Z_idx] = T_bits[i];
          Z_idx++;
        }
        nonsingletons_idx++;
      }
    }
    std::cout << "Y_str_idx: " << Y_str_idx << " Z_idx: " << Z_idx
              << " nonsingletons_idx: " << nonsingletons_idx << '\n';
    concat = Base4RankVectorTransposed<4>(Y_str);
    nonsingleton_sets = Pred8v2(nonsingleton_sets_vec);
    /* empty_sets = Pred8v2(empty_sets_plain); */
    empty_sets = bitvector_t(empty_sets_plain);
    sdsl::util::init_support(empty_sets_rs, &empty_sets);
    Z_A = bitvector_t(Z_A_plain);
    Z_C = bitvector_t(Z_C_plain);
    Z_G = bitvector_t(Z_G_plain);
    Z_T = bitvector_t(Z_T_plain);
    sdsl::util::init_support(Z_A_rs, &Z_A);
    sdsl::util::init_support(Z_C_rs, &Z_C);
    sdsl::util::init_support(Z_G_rs, &Z_G);
    sdsl::util::init_support(Z_T_rs, &Z_T);

    // For debugging: verify that ranks match
    /* rank_support_t A_bits_rs;
    rank_support_t C_bits_rs;
    rank_support_t G_bits_rs;
    rank_support_t T_bits_rs;

    sdsl::util::init_support(A_bits_rs, &(A_bits));
    sdsl::util::init_support(C_bits_rs, &(C_bits));
    sdsl::util::init_support(G_bits_rs, &(G_bits));
    sdsl::util::init_support(T_bits_rs, &(T_bits));

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
      }
      if (!(rC == ownC)) {
        std::cerr << "Rank mismatch at position " << i << " for C: " << rC
                  << " vs " << ownC << '\n';
      }
      if (!(rG == ownG)) {
        std::cerr << "Rank mismatch at position " << i << " for G: " << rG
                  << " vs " << ownG << '\n';
      }
      if (!(rT == ownT)) {
        std::cerr << "Rank mismatch at position " << i << " for T: " << rT
                  << " vs " << ownT << '\n';
      }
    } */
  }
};

}  // namespace sbwt