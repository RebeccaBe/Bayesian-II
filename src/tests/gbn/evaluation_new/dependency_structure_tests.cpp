#include "../../../../libs/catch/catch.hpp"

#include <fstream>
#include <iostream>
#include <chrono>

#include "../../../gbn/evaluation_new/wire_structure.h"
#include "../../../gbn/evaluation_new/evaluation.h"
#include "../../../gbn/general/gbn.h"
#include "../../../gbn/matrix/matrix.h"
#include "../../../gbn/matrix/matrix_io.h"
#include "../../../gbn/evaluation/wire_structure.h"
#include "../../../gbn/evaluation/evaluation.h"

#include "../../test_helpers.h"

/*TEST_CASE("wire structure tests") {
    auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f_matrices_3.gbn");
    auto wire_structure = improved::build_wire_structure(gbn);

    improved::print_wire_structure(std::cout, wire_structure, gbn);

}*/

TEST_CASE("comparing old and improved methods") {
    auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f_matrices_3.gbn");


    auto start_time = std::chrono::steady_clock::now();
    auto p_matrix_old = old::evaluate(gbn);
    auto end_time = std::chrono::steady_clock::now();
    double diff_milliseconds = std::chrono::duration<double, std::milli>(end_time-start_time).count();
    std::cout << "Old method:" << diff_milliseconds << std::endl;

    start_time = std::chrono::steady_clock::now();
    auto p_matrix_new = improved::evaluate(gbn);
    end_time = std::chrono::steady_clock::now();
    diff_milliseconds = std::chrono::duration<double, std::milli>(end_time-start_time).count();
    std::cout << "New method:" << diff_milliseconds << std::endl;


    REQUIRE(check_matrices_equal(*p_matrix_new,*p_matrix_old));
}

TEST_CASE("The evaluation of F-Matrices should result into a diagonal matrix.") {
    auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f_matrices.gbn");
    auto p_m = improved::evaluate(gbn);
    REQUIRE(p_m->type == DIAGONAL);
}

TEST_CASE("The evaluation of F-Matrices with swapped wires results into a non-diagonal matrix.") {
    auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f_matrices_swapped.gbn");
    auto p_m = old::evaluate(gbn); //can't be diganoal
}