#pragma once

#include "../general/gbn.h"

bool check_and_apply_F1(GBN& gbn, Vertex v, std::string& op);
bool check_and_apply_F2(GBN& gbn, Vertex v, std::string& op);
bool check_and_apply_CoUnit(GBN& gbn, Vertex v, std::string& op);
bool check_and_apply_F3(GBN& gbn, Vertex v, std::string& op);
bool check_and_apply_F4(GBN& gbn, Vertex v, std::string& op);
bool check_and_apply_F5(GBN& gbn, Vertex v, std::string& op);

bool split_vertex_if_multiple_outputs(GBN& gbn, Vertex v, std::string& op);
bool simplify_matrix_for_duplicate_inputs(GBN& gbn, Vertex v, std::string& op);
bool merge_F_matrices_to_diagonal_matrix(GBN& gbn, Vertex v, std::string& op);
bool reduce_diagonal_matrix(GBN& gbn, Vertex v, std::string& op);

std::set<Vertex> find_recursively_out(GBN gbn, Vertex v);
std::set<Vertex> find_recursively_in(GBN gbn, Vertex v);

