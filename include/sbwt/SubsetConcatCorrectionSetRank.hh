#pragma once

#include <map>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include <vector>

#include "Base4RankVectorTransposed.hh"
#include "Pred8v2.hh"
#include "Pred8v3.hh"
#include "Pred8vPino.hh"
#include "globals.hh"

namespace sbwt {

using namespace std;

template <typename bitvector_t, typename rank_support_t, typename pred8_cs_t>
class SubsetConcatCorrectionSetRank {
  Base4RankVectorTransposed<4>
      concat;          // The concatenated characters of all subsets
  Pred8vPino lengths;  // the positions in concat of sets larger than 1

  // bitvector_t T_corr_set;  // The universe of positions of T and $
  // rank_support_t T_corr_set_rs;

  pred8_cs_t T_corr_set;  // The universe of positions of T and $

 public:
  // Count of character c in subsets up to pos, not including pos
  int64_t rank(int64_t pos, char c) const {
    /* std::cout << "Rank called for pos " << pos << " and char " << (int)c
              << '\n'; */
    uint64_t pos_in_concat = lengths.select0(pos);
    /* uint64_t pos_in_concat_ = lengths.select0_Select(pos);
    if (pos_in_concat != pos_in_concat_) {
      cout << "Shit hit fan!" << endl;
      cout << " got " << pos_in_concat_ << " should got " << pos_in_concat <<
    endl; exit(1);
    } */
    /* if (pos < 10) {
      cout << "In rank searching for position pos = " << pos << " char " << c
           << " which corresponds to pos_in_concat = " << pos_in_concat << '\n';
    } */
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
    int64_t result = concat.rank(pos_in_concat, c_coded);
    if (c_coded == 3) {
      result -= T_corr_set.rank(result);
    }

    return result;
  }

  int64_t rank_by_charidx(int64_t pos, int64_t char_idx) {
    uint64_t pos_in_concat = lengths.select0(pos);

    uint64_t result = concat.rank(pos_in_concat, char_idx);
    if (char_idx == 3) {
      result -= T_corr_set.rank(result);
    }

    return result;
  }

  bool contains(int64_t pos, char c) const {
    // TODO: faster
    int64_t r1 = this->rank(pos, c);
    int64_t r2 = this->rank(pos + 1, c);
    return r1 != r2;
  }

  SubsetConcatCorrectionSetRank() {}

  SubsetConcatCorrectionSetRank(const sdsl::bit_vector& A_bits,
                                const sdsl::bit_vector& C_bits,
                                const sdsl::bit_vector& G_bits,
                                const sdsl::bit_vector& T_bits) {
    assert(A_bits.size() == C_bits.size() && C_bits.size() == G_bits.size() &&
           G_bits.size() == T_bits.size());

    int64_t n = A_bits.size();
    int64_t n_ones = 0;       // Number of the sets larger than 2 in SBWT
    int64_t n_empt = 0;       // Number of empty sets in SBWT
    int64_t sum_of_lens = 0;  // Number of symbols in SBWT
    int64_t T_universe_size = 0;
    for (int64_t i = 0; i < n; i++) {
      sum_of_lens += A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i];
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 0) {
        n_empt++;
      }
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] > 1)
        n_ones += (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] - 1);
      if (T_bits[i]) T_universe_size++;
    }
    T_universe_size += n_empt;
    std::vector<uint64_t> lens_positions;
    std::vector<uint64_t> T_corr_set_plain;
    // sdsl::bit_vector T_corr_plain(T_universe_size);
    std::string Y_str(sum_of_lens + n_empt, '\0');

    cout << "Constructing SubsetConcatCorrectionSetRank with n = " << n
         << " sum_of_lens = " << sum_of_lens << " n_ones = " << n_ones
         << " n_empt = " << n_empt << " T_universe_size = " << T_universe_size
         << '\n';
    int64_t Y_str_idx = 0, lens_idx = 0;
    uint64_t len_counter = 0;
    uint64_t t_counter = 0;
    for (int64_t i = 0; i < n; i++) {
      if (A_bits[i] == 1) Y_str[Y_str_idx++] = 0;
      if (C_bits[i] == 1) {
        Y_str[Y_str_idx++] = 1;
      }
      if (G_bits[i] == 1) Y_str[Y_str_idx++] = 2;
      if (T_bits[i] == 1) {
        Y_str[Y_str_idx++] = 3;
        // T_corr_plain[t_counter] = 0;
        t_counter++;
      }
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] == 0) {
        Y_str[Y_str_idx++] = 3;
        T_corr_set_plain.push_back(t_counter);
        t_counter++;
        // T_corr_plain[t_counter] = 1;
      }
      if (A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] > 1) {
        uint32_t set_len = A_bits[i] + C_bits[i] + G_bits[i] + T_bits[i] - 1;
        while (set_len) {
          lens_positions.push_back(len_counter);
          len_counter++;
          set_len--;
        }
      }
      len_counter++;
    }
    assert(lens_positions.size() == n_ones);
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

    concat = Base4RankVectorTransposed<4>(Y_str);
    lengths = Pred8vPino(lens_positions);
    // T_corr_set = bitvector_t(T_corr_plain);

    // sdsl::util::init_support(T_corr_set_rs, &T_corr_set);
    T_corr_set = pred8_cs_t(T_corr_set_plain);

    // For debugging
    /* rank_support_t A_bits_rs;
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
        std::cout << "subset " << i << '\n';
        std::cout << "set = " << A_bits[i] << C_bits[i] << G_bits[i]
                  << T_bits[i] << '\n';

        auto p = lengths.select0(i);
        std::cout << "select0 = " << p << '\n';

        std::cout << "Y[p-2..p+2]:\n";
        for (int64_t k = p - 2; k <= p + 2; ++k)
          std::cout << k << " -> " << (int)Y_str[k] << '\n';
      }
      if (!(rC == ownC)) {
        std::cerr << "Rank mismatch at position " << i << " for C: " << rC
                  << " vs " << ownC << '\n';
        wrongs++;
        std::cout << "subset " << i << '\n';
        std::cout << "set = " << A_bits[i] << C_bits[i] << G_bits[i]
                  << T_bits[i] << '\n';
        auto p_1 = lengths.select0(i - 1);
        std::cout << "select0(i-1) = " << p_1 << '\n';
        std::cout << "rank(select0(i-1)) = " << this->rank(i - 1, 'C') << '\n';
        auto p = lengths.select0(i);
        std::cout << "select0(i) = " << p << '\n';

        std::cout << "Y[p-2..p+2]:\n";
        for (int64_t k = p - 2; k <= p + 2; ++k)
          std::cout << k << " -> " << (int)Y_str[k] << '\n';
      }
      if (!(rG == ownG)) {
        std::cerr << "Rank mismatch at position " << i << " for G: " << rG
                  << " vs " << ownG << '\n';
        wrongs++;
        std::cout << "subset " << i << '\n';
        std::cout << "set =  " << A_bits[i] << C_bits[i] << G_bits[i]
                  << T_bits[i] << '\n';

        auto p = lengths.select0(i);
        std::cout << "select0 = " << p << '\n';

        std::cout << "Y[p-2..p+2]:\n";
        for (int64_t k = p - 2; k <= p + 2; ++k)
          std::cout << k << " -> " << (int)Y_str[k] << '\n';
      }
      if (!(rT == ownT)) {
        std::cerr << "Rank mismatch at position " << i << " for T: " << rT
                  << " vs " << ownT << '\n';
        wrongs++;
        std::cout << "subset " << i << '\n';
        std::cout << "set = " << A_bits[i] << C_bits[i] << G_bits[i]
                  << T_bits[i] << '\n';

        auto p = lengths.select0(i);
        std::cout << "select0 = " << p << '\n';

        std::cout << "Y[p-2..p+2]:\n";
        for (int64_t k = p - 2; k <= p + 2; ++k)
          std::cout << k << " -> " << (int)Y_str[k] << '\n';
      }
      if (wrongs > 20) {
        std::cerr << "More than 20 wrong ranks, stopping checking.\n";
        exit(1);
        break;
      }
    } */
  }

  int64_t serialize(ostream& os) const {
    int64_t written = 0;
    written += concat.serialize(os);
    written += lengths.serialize(os);
    int64_t tmp_written_ = T_corr_set.serialize(os);
    written += tmp_written_;
    // int64_t tmp_written = T_corr_set_rs.serialize(os);
    // written += tmp_written;

    std::cout << "Serialized concat " << concat.size_in_bytes() << " bytes\n";
    std::cout << "Serialized the lengths " << lengths.sizeInBytes()
              << " bytes\n";
    std::cout << "Serialized correction set of T " << tmp_written_
              << " bytes\n";
    /* std::cout << "Serialized rank support structure for correction set of T "
              << tmp_written << " bytes\n"; */
    return written;
  }

  void load(istream& is) {
    concat.load(is);
    lengths.load(is);
    T_corr_set.load(is);
    /* if (std::is_same<sdsl::rank_support_v5<>, rank_support_t>::value) {
      // Special handling needed for rank_support_v5 because of a design flaw in
      // sdsl
      T_corr_set_rs.load(is, &T_corr_set);
    } else {
      T_corr_set_rs.load(is);
      T_corr_set_rs.set_vector(&T_corr_set);
    } */

    cout << "Loaded SubsetConcatCorrectionSetRank with concat size: "
         << concat.size_in_bytes() << " lengths size: " << lengths.sizeInBytes()
         << " T corrections set size: " << T_corr_set.sizeInBytes() << '\n';
    cout << "For T correction set using " << typeid(T_corr_set).name() << '\n';
    cout << "For lengths using " << typeid(lengths).name() << '\n';
  }

  SubsetConcatCorrectionSetRank(const SubsetConcatCorrectionSetRank& other) {
    assert(&other != this);  // What on earth are you trying to do?
    operator=(other);
  }

  SubsetConcatCorrectionSetRank& operator=(
      const SubsetConcatCorrectionSetRank& other) {
    if (&other != this) {
      this->concat = other.concat;
      this->lengths = other.lengths;
      this->T_corr_set = other.T_corr_set;
      // this->T_corr_set_rs = other.T_corr_set_rs;
      // this->T_corr_set_rs.set_vector(&this->T_corr_set);
      return *this;
    } else
      return *this;  // Assignment to self -> do nothing.
  }
};

}  // namespace sbwt