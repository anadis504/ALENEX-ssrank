#pragma once

#include <map>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include <vector>

#include "Base4RankVector.hh"
#include "Base4RankVectorTransposed.hh"
#include "Base4RankVectorWordPacked.hh"
#include "Base4SummingRankVectorTransposed.hh"
#include "Pred8v2.hh"
#include "globals.hh"

namespace sbwt {

using namespace std;

template <typename bitvector_t, typename rank_support_t>
class SubsetNewConcatRank {
  Pred8v2 nonsingleton_sets;
  Base4SummingRankVectorTransposed<4> nonsingleton_lens;
  Base4RankVectorTransposed<4>
      concat;  // The concatenated characters of all subsets

 public:
  // Count of character c in subsets up to pos, not including pos
  int64_t rank(int64_t pos, char c) const {
    /* std::cout << "Rank called for pos " << pos << " and char " << (int)c
              << '\n'; */
    uint64_t nnz = nonsingleton_sets.rank(pos);
    uint64_t singleton_pos = pos - nnz;
    int64_t lens_sum = nonsingleton_lens.sum_of_ranks(nnz);
    // if (lens_sum)
    /* std::cout << "Pos: " << pos << " non-singletons: " << nnz
              << " Rank: " << rank << " Lens sum: " << lens_sum << '\n'; */
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
    int64_t result = concat.rank(lens_sum + singleton_pos, c_coded);

    return result;
  }

  int64_t rank_by_charidx(int64_t pos, int64_t char_idx) const {
    uint64_t nnz = nonsingleton_sets.rank(pos);
    uint64_t singleton_pos = pos - nnz;
    int64_t lens_sum = nonsingleton_lens.sum_of_ranks(nnz);

    int64_t result = concat.rank(lens_sum + singleton_pos, char_idx);

    return result;
  }

  bool contains(int64_t pos, char c) const {
    // TODO: faster
    int64_t r1 = this->rank(pos, c);
    int64_t r2 = this->rank(pos + 1, c);
    return r1 != r2;
  }

  SubsetNewConcatRank() {}

  SubsetNewConcatRank(const sdsl::bit_vector& A_bits,
                      const sdsl::bit_vector& C_bits,
                      const sdsl::bit_vector& G_bits,
                      const sdsl::bit_vector& T_bits) {
    assert(A_bits.size() == C_bits.size() && C_bits.size() == G_bits.size() &&
           G_bits.size() == T_bits.size());

    int64_t n = A_bits.size();
    int64_t n_b = 0;  // Number of branching nodes plus the nodes that do not
                      // have outedges
    int64_t n_u = 0;  // Number of nodes with exactly one outgoing edge
    int64_t sum_of_lens = 0;
    for (int64_t i = 0; i < n; i++) {
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 1)
        n_u++;
      else {
        n_b++;
        sum_of_lens += A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i];
      }
    }
    /* std::cout << "n: " << n << " n_b: " << n_b << " n_u: " << n_u
              << " sum_of_lens: " << sum_of_lens
              << " n_u+sum_of_lens: " << n_u + sum_of_lens << '\n'; */
    /* sdsl::bit_vector X_plain(n, 0); */
    std::vector<uint64_t> nonsingleton_sets_indices;
    std::string Y_str(sum_of_lens + n_u, '\0');
    std::string lens_vec(n_b, '\0');
    int64_t Y_str_idx = 0, lens_idx = 0;
    uint64_t lens[4] = {0, 0, 0, 0};
    for (int64_t i = 0; i < n; i++) {
      if (A_bits[i] == 1) Y_str[Y_str_idx++] = 0;
      if (C_bits[i] == 1) Y_str[Y_str_idx++] = 1;
      if (G_bits[i] == 1) Y_str[Y_str_idx++] = 2;
      if (T_bits[i] == 1) Y_str[Y_str_idx++] = 3;
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] != 1) {
        // One outgoing label
        /* X_plain[i] = 0;
      } else { */
        // 0 or > 1 outgoing labels
        /* X_plain[i] = 1; */
        nonsingleton_sets_indices.push_back(i);
        int64_t set_len = A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i];
        uint8_t encoded_len = (uint8_t)set_len - (uint8_t)(bool)(set_len > 0);
        lens[encoded_len]++;
        lens_vec[lens_idx++] = encoded_len;  // store length minus one
        /* if (encoded_len)
        std::cout << "Length for set " << i << " (" << lens_idx - 1 << ") : " <<
        set_len << " Encoded length for position " << i << ": "
                   << (int)encoded_len << '\n'; */
      }
    }
    /* for (int i = 0, j = 0; i < 20; i++) {
      std::cout << " len at pos " << i << ": " << (int)lens_vec[i];
      for (int k = 0; k < (int)lens_vec[i] + (bool)(int)lens_vec[i]; k++)
        std::cout << (int)Y_str[j++];
      std::cout << '\n';
    } */
    /* std::cout << "Constructed X of length " << X_plain.size()
              << " and Y of length " << Y_str.size()
              << " Y_str_idx: " << Y_str_idx << " lens_idx: " << lens_idx
              << '\n'; */
    /* nonsingleton_sets = bitvector_t(X_plain);
    sdsl::util::init_support(nonsingleton_sets_rs, &nonsingleton_sets); */
    nonsingleton_sets = Pred8v2(nonsingleton_sets_indices);
    nonsingleton_lens = Base4SummingRankVectorTransposed<4>(lens_vec);
    concat = Base4RankVectorTransposed<4>(Y_str);

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

  int64_t serialize(ostream& os) const {
    int64_t written = 0;
    written += concat.serialize(os);
    std::cout << "Serialized concat " << written << " bytes\n";
    written += nonsingleton_sets.serialize(os);
    std::cout << "Serialized nonsingleton_sets "
              << nonsingleton_sets.sizeInBytes() << " bytes\n";
    /* written += nonsingleton_sets_rs.serialize(os);
    std::cout << "Serialized nonsingleton_sets_rs " << written << " bytes\n"; */
    written += nonsingleton_lens.serialize(os);
    std::cout << "Serialized nonsingleton_lens "
              << nonsingleton_lens.size_in_bytes() << " bytes\n";
    return written;
  }

  void load(istream& is) {
    concat.load(is);
    nonsingleton_sets.load(is);

    /* if (std::is_same<sdsl::rank_support_v5<>, rank_support_t>::value) {
      nonsingleton_sets_rs.load(is, &nonsingleton_sets);
    } else {
      nonsingleton_sets_rs.load(is);
      nonsingleton_sets_rs.set_vector(&nonsingleton_sets);
    } */
    nonsingleton_lens.load(is);
    cout << "Loaded SubsetNewConcatRank with concat size: "
         << concat.size_in_bytes()
         << " nonsingleton_sets size: " << nonsingleton_sets.sizeInBytes()
         << " nonsingleton_lens size: " << nonsingleton_lens.size_in_bytes()
         << '\n';
    cout << "For nonsingelton_lens using " << typeid(nonsingleton_lens).name()
         << '\n';
    cout << "For nonsingelton_sets using " << typeid(nonsingleton_sets).name()
         << '\n';
  }

  SubsetNewConcatRank(const SubsetNewConcatRank& other) {
    assert(&other != this);  // What on earth are you trying to do?
    operator=(other);
  }

  SubsetNewConcatRank& operator=(const SubsetNewConcatRank& other) {
    if (&other != this) {
      this->concat = other.concat;
      this->nonsingleton_sets = other.nonsingleton_sets;
      /* this->nonsingleton_sets_rs = other.nonsingleton_sets_rs;
      this->nonsingleton_sets_rs.set_vector(&this->nonsingleton_sets); */
      this->nonsingleton_lens = other.nonsingleton_lens;
      return *this;
    } else
      return *this;  // Assignment to self -> do nothing.
  }
};

}  // namespace sbwt