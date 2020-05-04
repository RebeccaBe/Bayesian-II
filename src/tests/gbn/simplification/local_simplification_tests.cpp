#include "../../../../libs/catch/catch.hpp"

#include "../../../gbn/general/subgbn.h"
#include "../../../gbn/general/check.h"
#include "../../../gbn/simplification/simplification.h"
#include "../../../gbn/simplification/local_simplification.h"
#include "../../../gbn/evaluation/evaluation.h"
#include "../../../gbn/matrix/matrix_io.h"
#include "../../test_helpers.h"

#include <chrono>

TEST_CASE("OneB replacement (F1) should work correctly.")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "multiple_1_bs.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { simplification(gbn); return gbn; });
}

TEST_CASE("(F2) should work correctly.")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f2_simplification.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { simplification(gbn);simplification(gbn);simplification(gbn); return gbn; });
}

TEST_CASE("(F3) should work correctly (f3_simplification_1.gbn).")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f3_simplification_1.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; check_and_apply_F3(gbn,5,s); return gbn; });
}

TEST_CASE("(F3) should work correctly (f3_simplification_2.gbn).")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f3_simplification_2.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; check_and_apply_F3(gbn,5,s); return gbn; });
}

TEST_CASE("(F3) should work correctly (f3_simplification_3.gbn).")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f3_simplification_3.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; check_and_apply_F3(gbn,6,s); return gbn; });
}

TEST_CASE("(F3) should work correctly (f3_simplification_4.gbn).")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f3_simplification_4.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; check_and_apply_F3(gbn,2,s); return gbn; });
}

TEST_CASE("(F4) should work correctly (f4_simplification_1.gbn).")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f4_simplification_1.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; check_and_apply_F4(gbn,0,s); return gbn; }, [](GBN gbn_before, GBN gbn_after) -> void {
		REQUIRE(boost::num_edges(gbn_after.graph) == boost::num_edges(gbn_before.graph)-1);	
	});
}

TEST_CASE("(F4) should work correctly (f4_simplification_2.gbn).")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f4_simplification_2.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; check_and_apply_F4(gbn,0,s); return gbn; }, [](GBN gbn_before, GBN gbn_after) -> void {
		REQUIRE(boost::num_edges(gbn_after.graph) == boost::num_edges(gbn_before.graph)-1);	
	});
}

// Diagonal Matrices changed
/*TEST_CASE("(F4) should work correctly with corresponding diagonal matrix (f4_simplification_1D.gbn).")
{
    auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f4_simplification_1D.gbn");
    check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; check_and_apply_F4(gbn,0,s); return gbn; }, [](GBN gbn_before, GBN gbn_after) -> void {
        REQUIRE(boost::num_edges(gbn_after.graph) == boost::num_edges(gbn_before.graph)-1);
    });
}

TEST_CASE("(F4) should work correctly with corresponding diagonal matrix (f4_simplification_2D.gbn).")
{
    auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f4_simplification_2D.gbn");
    check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; check_and_apply_F4(gbn,0,s); return gbn; }, [](GBN gbn_before, GBN gbn_after) -> void {
        REQUIRE(boost::num_edges(gbn_after.graph) == boost::num_edges(gbn_before.graph)-1);
    });
}*/

TEST_CASE("(F5) should work correctly (f5_simplification_1.gbn).")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f5_simplification_1.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; check_and_apply_F5(gbn,0,s); return gbn; }, [](GBN gbn_before, GBN gbn_after) -> void {
		REQUIRE(boost::num_edges(gbn_before.graph) == 8);	
		REQUIRE(boost::num_edges(gbn_after.graph) == 5);	
	});
}

/*TEST_CASE("(F5) should work correctly with corresponding diagonal matrix (f5_simplification_1D.gbn).")
{
    auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f5_simplification_1D.gbn");
    check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; check_and_apply_F5(gbn,0,s); return gbn; }, [](GBN gbn_before, GBN gbn_after) -> void {
        REQUIRE(boost::num_edges(gbn_before.graph) == 8);
        REQUIRE(boost::num_edges(gbn_after.graph) == 5);
    });
}*/

TEST_CASE("duplicate_inputs.gbn: (simplify_matrix_for_duplicate_inputs)")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "duplicate_inputs.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; simplify_matrix_for_duplicate_inputs(gbn,0,s); return gbn; });
}

TEST_CASE("duplicate_inputs2.gbn: (simplify_matrix_for_duplicate_inputs)")
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "duplicate_inputs2.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; simplify_matrix_for_duplicate_inputs(gbn,0,s); return gbn; });
}

TEST_CASE("F matrices merging to diagonal matrix should work correctly (f_matrices).") {

	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f_matrices.gbn");
	auto p_m_before = evaluate(gbn);

    check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; merge_F_matrices_to_diagonal_matrix(gbn,0,s); return gbn; }, [](GBN gbn_before, GBN gbn_after) -> void {
        REQUIRE(boost::num_edges(gbn_before.graph) == 17);
        REQUIRE(boost::num_edges(gbn_after.graph) == 10);
        REQUIRE(matrix(inside_vertices(gbn_after).at(0), gbn_after.graph)->type == DIAGONAL);
        REQUIRE(inside_vertices(gbn_after).size() == 1);
    });
}

TEST_CASE("independent_vars.gbn: (reducing a diagonal matrix with independency)") {

	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "independent_vars.gbn");
	check_evaluates_equal_after_operation(gbn, [](GBN gbn) -> GBN { std::string s; reduce_diagonal_matrix(gbn,0,s); return gbn; }, [](GBN gbn_before, GBN gbn_after) -> void {
		REQUIRE(boost::num_edges(gbn_before.graph) == 7);
		REQUIRE(boost::num_edges(gbn_after.graph) == 6);
		REQUIRE(matrix(inside_vertices(gbn_after).at(0), gbn_after.graph)->n == 2);
	});
}