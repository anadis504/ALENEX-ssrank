#pragma once

#include <vector>
#include <string>

std::vector<std::string> get_available_variants();

int build_main(int argc, char** argv);
int search_main(int argc, char** argv);
int build_from_plain_main(int argc, char** argv);
int ascii_export_main(int argc, char** argv);
int s_rank_queries_main(int argc, char** argv);
int sample_kmers_from_index_main(int argc, char** argv);
