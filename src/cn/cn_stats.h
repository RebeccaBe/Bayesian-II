#include "cn.h"


std::size_t get_n_valid_transitions(CN cn);
int get_n_transitions(CN cn);
int get_n_markings(CN cn);
std::vector<std::size_t> get_distribution_pre(CN cn);
std::vector<std::size_t> get_distribution_post(CN cn);
double get_average_pre(CN cn);
double get_average_post(CN cn);
std::vector<std::size_t> get_valid_transitions(CN cn);