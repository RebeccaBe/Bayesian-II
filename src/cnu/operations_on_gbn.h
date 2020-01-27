#pragma once

#include "../gbn/general/gbn.h"
#include "../helpers.hpp"
#include "../gbn/matrix/matrix_io.h"

void set_op(const std::vector<std::size_t> places, bool b, GBN& gbn);
void assert_op(const std::vector<std::size_t> places, bool b, GBN& gbn);
void nassert_op(const std::vector<std::size_t> places, bool b, GBN& gbn);
void setp_op(const std::vector<std::size_t> places, double p, GBN& gbn);
void successp_op(const std::vector<std::vector<std::size_t>> pre_places,
                 const std::vector<std::vector<std::size_t>> post_places,
                 const std::vector<double> probabilities, GBN& gbn);
void successStoch_op(const std::vector<std::vector<std::size_t>> pre_places,
                 const std::vector<std::vector<std::size_t>> post_places,
                 const std::vector<double> probabilities, GBN& gbn);

bool validate_transition(std::vector<std::size_t> places, bool condition, GBN& gbn);
double normalization_factor(BitVec assignmentX, std::map<std::size_t, std::size_t> mapping_place_key,
        const std::vector<std::vector<std::size_t>> pre_places,
        const std::vector<double> probabilities); // post-places not needed due to no post-condition
void normalize_matrix_rows(Matrix& matrix);
