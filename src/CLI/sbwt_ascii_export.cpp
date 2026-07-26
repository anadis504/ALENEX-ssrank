#include "SBWT.hh"
#include "SeqIO/SeqIO.hh"
#include "SubsetMatrixRank.hh"
#include "commands.hh"
#include "cxxopts.hpp"
#include "globals.hh"
#include "throwing_streams.hh"
#include "variants.hh"

using namespace std;
using namespace sbwt;

template <typename sbwt_t>
void export_sbwt_variant(const sbwt_t& sbwt, seq_io::Buffered_ofstream<>& out) {
  sbwt.ascii_export_metadata(out);
  out.flush();
  sbwt.ascii_export_sets(out);
  out.flush();
}

int ascii_export_main(int argc, char** argv) {
  cxxopts::Options options(
      argv[0],
      "\nExport SBWT in an ASCII format. First comes three lines of "
      "metadata:\n\nversion: <version>\nk: <k>\nnumber_of_sets: <number of "
      "sets>\nnumber_of_kmers: <number of k-mers>\n\nThis is followed by a "
      "line of the form \"sbwt: <sbwt>\". The format of <sbwt> is as follows: "
      "Each set is written as a string of characters. A non-empty set is a "
      "string of DNA characters (ACGT), such that the last character of the "
      "set is written in lower case. Empty sets are represented as a single "
      "'$'. The representations of the sets are concatenated together. For "
      "example, the sequence {A,C}, {A,T}, {}, {C}, {A,G,T} is represented as "
      "AcAt$cAGt.");

  options.add_options()("o,out-file", "Output filename.",
                        cxxopts::value<string>())(
      "i,index-file", "Index input file.", cxxopts::value<string>())(
      "h,help", "Print usage");

  int64_t old_argc = argc;  // Must store this because the parser modifies it
  auto opts = options.parse(argc, argv);

  if (old_argc == 1 || opts.count("help")) {
    std::cerr << options.help() << std::endl;
    exit(1);
  }

  string indexfile = opts["index-file"].as<string>();
  string outfile = opts["out-file"].as<string>();

  vector<string> variants = get_available_variants();

  throwing_ifstream in(indexfile, ios::binary);
  string variant = load_string(in.stream);  // read variant type
  if (std::find(variants.begin(), variants.end(), variant) == variants.end()) {
    cerr << "Error loading index from file: unrecognized variant specified in "
            "the file"
         << endl;
    return 1;
  }

  write_log("Loading and exporting the index variant " + variant,
            LogLevel::MAJOR);
  seq_io::Buffered_ofstream<> out(outfile);

  if (variant == "plain-matrix") {
    plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }
  if (variant == "mef-matrix") {
    cerr << "Error: Index export does not work for mef-matrix because mef does "
            "not implement access to the sets"
         << endl;
    return 1;
  }
  if (variant == "plain-split") {
    plain_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }
  if (variant == "mef-split") {
    cerr << "Error: Index export does not work for mef-split because mef does "
            "not implement access to the sets"
         << endl;
    return 1;
  }
  if (variant == "correction-sets") {
    correction_sets_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }
  if (variant == "npred8-split-w-packed") {
    new_split_w_packed_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }
  if (variant == "blocked-correction-sets") {
    blocked_correction_sets_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }
  if (variant == "fixed-block-correction-setsA") {
    fixed_block_correction_sets1_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }
  if (variant == "fixed-block-correction-setsB") {
    fixed_block_correction_sets2_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }

  if (variant == "fixed-block-correction-setsC") {
    fixed_block_correction_sets3_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }
  if (variant == "blocked8-split") {
    blocked8_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }
  if (variant == "blocked9-split") {
    blocked9_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }
  if (variant == "pred8-WT-split") {
    pred8_split_sbwt_t sbwt;
    sbwt.load(in.stream);
    export_sbwt_variant(sbwt, out);
  }

  return 0;
}
