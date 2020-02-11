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
void normalize_matrix_rows(Matrix& matrix);
void normalize_matrix_cols(Matrix& matrix);

