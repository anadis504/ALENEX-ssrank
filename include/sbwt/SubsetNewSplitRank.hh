#pragma once

#include <map>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <vector>

#include "Base4RankVector.hh"
#include "Base4RankVectorTransposed.hh"
#include "Base4RankVectorWordPacked.hh"
#include "Pred16.hh"
#include "Pred8v2.hh"
#include "SBWT.hh"
#include "globals.hh"

namespace sbwt {

/* template <typename X_bitvector_t, typename X_bitvector_rank_t,
          typename Z_bitvector_t, typename Z_rank_support_t> */
template <typename X_pred8_t, typename Y_quadrankvector_t,
          typename Z_bitvector_t, typename Z_rank_support_t>
class SubsetNewSplitRank {
 public:
  X_pred8_t X;

  Y_quadrankvector_t Y;
  // Base4RankVectorTransposed<4> Y;

  Z_bitvector_t
      Z_A;  // The row of 'A' in the Z matrix, in columns with != 1 one-bit
  Z_bitvector_t
      Z_C;  // The row of 'C' in the Z matrix, in columns with != 1 one-bit
  Z_bitvector_t
      Z_G;  // The row of 'G' in the Z matrix, in columns with != 1 one-bit
  Z_bitvector_t
      Z_T;  // The row of 'T' in the Z matrix, in columns with != 1 one-bit

  Z_rank_support_t* rank_supports[4];
  Z_rank_support_t Z_A_rs;
  Z_rank_support_t Z_C_rs;
  Z_rank_support_t Z_G_rs;
  Z_rank_support_t Z_T_rs;
  // Returns the number of bytes written
  int64_t serialize(ostream& os) const {
    int64_t written = 0;
    int64_t written_tmp = 0;
    written += X.serialize(os);
    cerr << "X serialized " << written << " bytes\n";
    written_tmp = Y.serialize(os);
    cerr << "Y serialized " << written_tmp << " bytes\n";
    written += written_tmp;

    written_tmp = Z_A.serialize(os);
    cerr << "Z_A serialized " << written_tmp << " bytes\n";
    written += written_tmp;
    written += Z_C.serialize(os);
    written += Z_G.serialize(os);
    written += Z_T.serialize(os);

    // written_tmp = X_rs.serialize(os);
    // written += written_tmp;
    // cerr << "X_rs serialized " << written_tmp << " bytes\n";
    written_tmp  = Z_A_rs.serialize(os);
    cerr << "Z_A_rs serialized " << written_tmp << " bytes\n";
    written += written_tmp;
    written += Z_C_rs.serialize(os);
    written += Z_G_rs.serialize(os);
    written += Z_T_rs.serialize(os);

    return written;
  }

  void load(istream& is) {
    X.load(is);
    Y.load(is);

    Z_A.load(is);
    Z_C.load(is);
    Z_G.load(is);
    Z_T.load(is);

    if (std::is_same<sdsl::rank_support_v5<>, Z_rank_support_t>::value) {
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
    this->rank_supports[0] = &(this->Z_A_rs);
    this->rank_supports[1] = &(this->Z_C_rs);
    this->rank_supports[2] = &(this->Z_G_rs);
    this->rank_supports[3] = &(this->Z_T_rs);
    cerr << "uses " << typeid(Y).name() << "\n";
  }

  SubsetNewSplitRank() {}

  SubsetNewSplitRank(const sdsl::bit_vector& A_bits,
                     const sdsl::bit_vector& C_bits,
                     const sdsl::bit_vector& G_bits,
                     const sdsl::bit_vector& T_bits) {
    int64_t n = A_bits.size();
    int64_t n_b = 0;  // Number of branching nodes plus the nodes that do not
                      // have outedges
    int64_t n_u = 0;  // Number of nodes with exactly one outgoing edge

    std::vector<uint64_t> nonsingleton_sets_vec;
    for (int64_t i = 0; i < n; i++) {
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 1)
        n_u++;
      else {
        nonsingleton_sets_vec.push_back(i);
        n_b++;
      }
    }

    std::cerr << "n n_b n_u " << n << " " << n_b << " " << n_u << endl;

    sdsl::bit_vector X_plain(n);
    sdsl::bit_vector Z_A_plain(n_b);
    sdsl::bit_vector Z_C_plain(n_b);
    sdsl::bit_vector Z_G_plain(n_b);
    sdsl::bit_vector Z_T_plain(n_b);
    std::string Y_str(n_u, '\0');

    for (int64_t i = 0, Y_str_idx = 0, Z_idx = 0; i < n; i++) {
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 1) {
        // One outgoing label
        X_plain[i] = 0;
        if (A_bits[i] == 1) Y_str[Y_str_idx] = 0;  //'A';
        if (C_bits[i] == 1) Y_str[Y_str_idx] = 1;  //'C';
        if (G_bits[i] == 1) Y_str[Y_str_idx] = 2;  //'G';
        if (T_bits[i] == 1) Y_str[Y_str_idx] = 3;  //'T';
        Y_str_idx++;
      } else {
        // 0 or > 1 outgoing labels
        X_plain[i] = 1;
        Z_A_plain[Z_idx] = A_bits[i];
        Z_C_plain[Z_idx] = C_bits[i];
        Z_G_plain[Z_idx] = G_bits[i];
        Z_T_plain[Z_idx] = T_bits[i];
        Z_idx++;
      }
    }

    // X = Y_quadrankvector_type(X_plain);
    X = X_pred8_t(nonsingleton_sets_vec);

    // sdsl::construct_im(Y, Y_str.c_str(), 1); // 1: file format is a sequence,
    // not a serialized sdsl object
    Z_A = Z_bitvector_t(Z_A_plain);
    Z_C = Z_bitvector_t(Z_C_plain);
    Z_G = Z_bitvector_t(Z_G_plain);
    Z_T = Z_bitvector_t(Z_T_plain);

    // sdsl::util::init_support(X_rs, &X);
    sdsl::util::init_support(Z_A_rs, &Z_A);
    sdsl::util::init_support(Z_C_rs, &Z_C);
    sdsl::util::init_support(Z_G_rs, &Z_G);
    sdsl::util::init_support(Z_T_rs, &Z_T);

    this->rank_supports[0] = &(this->Z_A_rs);
    this->rank_supports[1] = &(this->Z_C_rs);
    this->rank_supports[2] = &(this->Z_G_rs);
    this->rank_supports[3] = &(this->Z_T_rs);

    Y = Y_quadrankvector_t(Y_str);

    // For debugging: verify that ranks match
    Z_rank_support_t A_bits_rs;
    Z_rank_support_t C_bits_rs;
    Z_rank_support_t G_bits_rs;
    Z_rank_support_t T_bits_rs;

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
    } 
  }

  SubsetNewSplitRank(const SubsetNewSplitRank& other) {
    assert(&other != this);  // What on earth are you trying to do?
    operator=(other);
  }

  SubsetNewSplitRank& operator=(const SubsetNewSplitRank& other) {
    if (&other != this) {
      this->X = other.X;
      this->Y = other.Y;

      this->Z_A = other.Z_A;
      this->Z_C = other.Z_C;
      this->Z_G = other.Z_G;
      this->Z_T = other.Z_T;

      // this->X_rs = other.X_rs;
      this->Z_A_rs = other.Z_A_rs;
      this->Z_C_rs = other.Z_C_rs;
      this->Z_G_rs = other.Z_G_rs;
      this->Z_T_rs = other.Z_T_rs;

      // this->X_rs.set_vector(&(this->X));
      this->Z_A_rs.set_vector(&(this->Z_A));
      this->Z_C_rs.set_vector(&(this->Z_C));
      this->Z_G_rs.set_vector(&(this->Z_G));
      this->Z_T_rs.set_vector(&(this->Z_T));

      this->rank_supports[0] = &(this->Z_A_rs);
      this->rank_supports[1] = &(this->Z_C_rs);
      this->rank_supports[2] = &(this->Z_G_rs);
      this->rank_supports[3] = &(this->Z_T_rs);

      return *this;
    } else
      return *this;  // Assignment to self -> do nothing.
  }

  int64_t rank(int64_t pos, char c) const {
    // std::cerr << "pos c: " << pos << ' ' << c << '\n';
    int64_t rank1 = X.rank(pos);
    int64_t rank0 = pos - rank1;
    uint8_t c_coded = DNA_to_char_idx(c);
    int64_t Y_count = Y.rank(rank0, c_coded);
    int64_t correction = 0;
    switch (c) {
      case 'A':

        correction = Z_A_rs.rank(rank1);
        break;
      case 'C':

        correction = Z_C_rs.rank(rank1);
        break;
      case 'G':

        correction = Z_G_rs.rank(rank1);
        break;
      case 'T':

        correction = Z_T_rs.rank(rank1);
        break;
      default:
        cerr << "Error: Rank called with non-ACGT character: " << c << endl;
        exit(1);
    }
    return Y_count + correction;
  }

  int64_t rank_by_charidx(int64_t pos, int64_t char_idx) const {
    int64_t rank1 = X.rank(pos);
    int64_t rank0 = pos - rank1;

    int64_t Y_count = Y.rank(rank0, char_idx_to_DNA(char_idx));
    int64_t correction = (*rank_supports[char_idx]).rank(rank1);

    return Y_count + correction;
  }

  bool contains(int64_t pos, char c) const {
    // TODO: faster
    int64_t r1 = this->rank(pos, c);
    int64_t r2 = this->rank(pos + 1, c);
    return r1 != r2;
  }
};

}  // namespace sbwt
