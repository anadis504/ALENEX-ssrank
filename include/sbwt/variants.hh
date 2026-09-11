#pragma once

#include <cstring>
#include <filesystem>
#include <string>

#include "Base4RankVector.hh"
#include "Base4RankVectorTransposed.hh"
#include "Base4RankVectorWordPacked.hh"
#include "BlockedCorrectionSetsConstant.hh"
#include "BlockedCorrectionSetsConstantSmaller.hh"
#include "BlockedCorrectionSetsConstantWordPacked6.hh"
#include "BlockedCorrectionSetsConstantWordPacked6smaller.hh"
#include "BlockedCorrectionSetsConstantWordPacked7.hh"
#include "BlockedCorrectionSetsConstantWordPacked7smaller.hh"
#include "MEF.hpp"
#include "Pred8vPinoLight.hh"
#include "SBWT.hh"
#include "SubsetBlockedCorrectionSets.hh"
#include "SubsetBlockedSplit8.hh"
#include "SubsetBlockedSplit9.hh"
#include "SubsetConcatCorrectionSetRank.hh"
#include "SubsetConcatSplitLenthsRank.hh"
#include "SubsetCorrectionSets.hh"
#include "SubsetFixedBlockCorrectionSets6.hh"
#include "SubsetFixedBlockCorrectionSets7.hh"
#include "SubsetMatrixRank.hh"
#include "SubsetNewConcatRank.hh"
#include "SubsetNewSplitRank.hh"
#include "SubsetSplitRank.hh"
#include "SubsetSplitRankPred8.hh"
#include "SubsetSplitSmallerSizeRank.hh"
#include "cxxopts.hpp"
#include "globals.hh"
#include "stdlib_printing.hh"

namespace sbwt {

// matrices
typedef SBWT<SubsetMatrixRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    plain_matrix_sbwt_t;

typedef SBWT<SubsetMatrixRank<mod_ef_vector<>, mod_ef_vector<>::rank_1_type>>
    mef_matrix_sbwt_t;  // Currently does not support extracting all k-mers
                        // because mod_ef_vector does not support access.

// splits
typedef SBWT<SubsetSplitRank<sdsl::bit_vector, sdsl::rank_support_v5<>,
                             sdsl::bit_vector, sdsl::rank_support_v5<>>>
    plain_split_sbwt_t;

typedef SBWT<SubsetSplitRank<mod_ef_vector<>, mod_ef_vector<>::rank_1_type,
                             sdsl::bit_vector, sdsl::rank_support_v5<>>>
    mef_split_sbwt_t;  // Currently does not support extracting all k-mers
                       // because mod_ef_vector does not support access.

typedef SBWT<SubsetSplitRank<sdsl::sd_vector<>, sdsl::sd_vector<>::rank_1_type,
                             sdsl::bit_vector, sdsl::rank_support_v5<>>>
    ef_split_sbwt_t;

// new splits
typedef SBWT<
    SubsetSplitRankPred8<Pred8v2, sdsl::bit_vector, sdsl::rank_support_v5<>>>
    pred8_split_sbwt_t;
typedef SBWT<SubsetSplitRankPred8<Pred8vPinoLight, sdsl::bit_vector,
                                  sdsl::rank_support_v5<>>>
    pred8_pino_split_sbwt_t;

typedef SBWT<SubsetNewSplitRank<Pred8v2, Base4RankVector<4>, sdsl::bit_vector,
                                sdsl::rank_support_v5<>>>
    new_split_packed_sbwt_t;

typedef SBWT<SubsetNewSplitRank<Pred8v2, Base4RankVectorWordPacked<4>,
                                sdsl::bit_vector, sdsl::rank_support_v5<>>>
    new_split_w_packed_sbwt_t;

typedef SBWT<SubsetNewSplitRank<Pred8v2, Base4RankVectorTransposed<4>,
                                sdsl::bit_vector, sdsl::rank_support_v5<>>>
    new_split_transposed_sbwt_t;

typedef SBWT<SubsetNewSplitRank<Pred8vPinoLight, Base4RankVectorTransposed<4>,
                                sdsl::bit_vector, sdsl::rank_support_v5<>>>
    new_split_pino_transposed_sbwt_t;

typedef SBWT<SubsetBlockedSplitRank8<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    blocked8_split_sbwt_t;
typedef SBWT<SubsetBlockedSplitRank9<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    blocked9_split_sbwt_t;

// Correction set
typedef SBWT<SubsetCorrectionSetsRank<Base4RankVectorTransposed<4>,
                                      sdsl::bit_vector, sdsl::rank_support_v5<>,
                                      Pred8vPinoLight>>
    correction_sets_sbwt_t;

typedef SBWT<
    SubsetBlockedCorrectionSetsRank<sdsl::bit_vector, sdsl::rank_support_v5<>>>
    blocked_correction_sets_sbwt_t;
typedef SBWT<SubsetFixedBlockCorrectionSetsRank7<
    FixedBlockedCorrectionSetsBase4Rank1<4>, sdsl::bit_vector,
    sdsl::rank_support_v5<>>>
    fixed_block_correction_sets1_sbwt_t;

typedef SBWT<SubsetFixedBlockCorrectionSetsRank7<
    FixedBlockedCorrectionSetsBase4Rank1_<4>, sdsl::bit_vector,
    sdsl::rank_support_v5<>>>
    fixed_block_correction_sets1_smaller_sbwt_t;

typedef SBWT<SubsetFixedBlockCorrectionSetsRank7<
    FixedBlockedCorrectionSetsBase4Rank2<4>, sdsl::bit_vector,
    sdsl::rank_support_v5<>>>
    fixed_block_correction_sets2_sbwt_t;

typedef SBWT<SubsetFixedBlockCorrectionSetsRank3<
    BlockedCorrectionSetsBase4Rank67<4>, sdsl::bit_vector,
    sdsl::rank_support_v5<>>>
    fixed_block_correction_sets3_sbwt_t;

}  // namespace sbwt
