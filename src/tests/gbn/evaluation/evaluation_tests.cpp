#include "../../../../libs/catch/catch.hpp"

#include <fstream>
#include <chrono>

#include "../../../gbn/evaluation/evaluation.h"
#include "../../../gbn/matrix/matrix_io.h"
#include "../../../gbn/general/check.h"
#include "../../../gbn/general/gbn_io.h"
#include "../../../gbn/general/subgbn.h"
#include "../../../gbn/modification/merging.h"
#include "../../../gbn/general/special_cases.h"
#include "../../../cn/randomized_generation.h"

#include "../../test_helpers.h"
#include "../../../runtime_tests/helpers/random_transition_helper.h"
#include "../../../cnu/fire_transition.h"
#include "../../../gbn/simplification/simplification.h"


TEST_CASE("Evaluation of single node should just give that matrix of node 1") 
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "one_node.gbn");

	auto sub_gbn = SubGBN::make_from_vertices(gbn, { 0 });
	auto p_m = evaluate(sub_gbn.gbn);

	REQUIRE(p_m->get(BitVec("0"),BitVec("0")) == Approx(0.333333333333));
	REQUIRE(p_m->get(BitVec("1"),BitVec("0")) == Approx(0.666666666666));
	REQUIRE(p_m->get(BitVec("0"),BitVec("1")) == Approx(0.75));
	REQUIRE(p_m->get(BitVec("1"),BitVec("1")) == Approx(0.25));
}

TEST_CASE("Evaluation of single node should just give that matrix of node 2") 
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "one_node_two_wires.gbn");

	std::vector<Vertex> vec({ 0 });

	auto sub_gbn = SubGBN::make_from_vertices(gbn, { 0 });
	auto wire_structure = build_wire_structure(sub_gbn.gbn);
	REQUIRE(wire_structure.wires.size() == 4);
	for(auto w : wire_structure.wires)
	{
		REQUIRE(w.inside_ports.size() == 1);
		REQUIRE(w.io_ports.size() == 1);
	}
	REQUIRE(sub_gbn.gbn.n == 2);
	REQUIRE(sub_gbn.gbn.m == 2);

	auto p_m = evaluate(gbn);

	REQUIRE(p_m->get(BitVec("00"),BitVec("00")) == Approx(0.25));
	REQUIRE(p_m->get(BitVec("00"),BitVec("01")) == Approx(0.333333333333));
	REQUIRE(p_m->get(BitVec("00"),BitVec("10")) == Approx(0.5));
	REQUIRE(p_m->get(BitVec("00"),BitVec("11")) == Approx(1.0));

	REQUIRE(p_m->get(BitVec("01"),BitVec("00")) == Approx(0.25));
	REQUIRE(p_m->get(BitVec("01"),BitVec("01")) == Approx(0.333333333333));
	REQUIRE(p_m->get(BitVec("01"),BitVec("10")) == Approx(0.5));
	REQUIRE(p_m->get(BitVec("01"),BitVec("11")) == Approx(0));

	REQUIRE(p_m->get(BitVec("10"),BitVec("00")) == Approx(0.25));
	REQUIRE(p_m->get(BitVec("10"),BitVec("01")) == Approx(0.333333333333));
	REQUIRE(p_m->get(BitVec("10"),BitVec("10")) == Approx(0));
	REQUIRE(p_m->get(BitVec("10"),BitVec("11")) == Approx(0));

	REQUIRE(p_m->get(BitVec("11"),BitVec("00")) == Approx(0.25));
	REQUIRE(p_m->get(BitVec("11"),BitVec("01")) == Approx(0));
	REQUIRE(p_m->get(BitVec("11"),BitVec("10")) == Approx(0));
	REQUIRE(p_m->get(BitVec("11"),BitVec("11")) == Approx(0));

}

TEST_CASE("For a single node gbn every wire should have exactly one v_input bitvec and one io bitvec") 
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "one_node_two_wires.gbn");

	auto sub_gbn = SubGBN::make_from_vertices(gbn, { 0 });
	auto wire_structure = build_wire_structure(sub_gbn.gbn);

	for(auto& w : wire_structure.wires)
	{
		REQUIRE(w.inside_ports.size() == 1);
		REQUIRE(w.io_ports.size() == 1);
	}

}

TEST_CASE("Evaluation for split.gbn should be correct") 
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "split.gbn");

	auto sub_gbn = SubGBN::make_from_vertices(gbn, { 0, 1, 2 });
	auto wire_structure = build_wire_structure(sub_gbn.gbn);

	auto p_m = evaluate(sub_gbn.gbn);
	
	// check dimensions
	REQUIRE(p_m->n == 1);
	REQUIRE(p_m->m == 2);

	// checkp_ values
	REQUIRE(p_m->get(BitVec("00"), BitVec("0")) == Approx(0.2037037037));
	REQUIRE(p_m->get(BitVec("01"), BitVec("0")) == Approx(0.2407407407));
	REQUIRE(p_m->get(BitVec("10"), BitVec("0")) == Approx(0.2407407407));
	REQUIRE(p_m->get(BitVec("11"), BitVec("0")) == Approx(0.3148148148));

	REQUIRE(p_m->get(BitVec("00"), BitVec("1")) == Approx(0.1805555556));
	REQUIRE(p_m->get(BitVec("01"), BitVec("1")) == Approx(0.2361111111));
	REQUIRE(p_m->get(BitVec("10"), BitVec("1")) == Approx(0.2361111111));
	REQUIRE(p_m->get(BitVec("11"), BitVec("1")) == Approx(0.3472222222));
}

TEST_CASE("Evaluation for four_nodes.gbn should be correct") 
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "four_nodes.gbn");

	auto sub_gbn = SubGBN::make_from_vertices(gbn, { 0, 1, 2, 3 });
	auto wire_structure = build_wire_structure(sub_gbn.gbn);
	auto p_m = evaluate(sub_gbn.gbn);
	
	// check dimensions
	REQUIRE(p_m->n == 2);
	REQUIRE(p_m->m == 1);

	// checkp_ values
	REQUIRE(p_m->get(BitVec("0"), BitVec("00")) == Approx(0.4074074074));
	REQUIRE(p_m->get(BitVec("0"), BitVec("01")) == Approx(0.4166666667));
	REQUIRE(p_m->get(BitVec("0"), BitVec("10")) == Approx(0.4305555556));
	REQUIRE(p_m->get(BitVec("0"), BitVec("11")) == Approx(0.4375000000));

	REQUIRE(p_m->get(BitVec("1"), BitVec("00")) == Approx(0.5925925926));
	REQUIRE(p_m->get(BitVec("1"), BitVec("01")) == Approx(0.5833333333));
	REQUIRE(p_m->get(BitVec("1"), BitVec("10")) == Approx(0.5694444444));
	REQUIRE(p_m->get(BitVec("1"), BitVec("11")) == Approx(0.5625000000));
}

TEST_CASE("Evaluation for seven_nodes.gbn should work.") 
{
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "seven_nodes.gbn");

	auto sub_gbn = SubGBN::make_from_vertices(gbn, { 0, 1, 2, 3, 4, 5, 6 });
	auto wire_structure = build_wire_structure(sub_gbn.gbn);
	auto wire_structure2 = build_wire_structure(gbn);
	auto p_m = evaluate(sub_gbn.gbn);

	REQUIRE(wire_structure.wires.size() == wire_structure2.wires.size());
	REQUIRE(wire_structure.wires.size() == 13);

	REQUIRE(p_m->n == 4);
	REQUIRE(p_m->m == 3);
	REQUIRE(is_stochastic(*p_m) == true);
}

TEST_CASE("Evaluating id_1.gbn should yield identity 1 -> 1") {
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "id_1.gbn");

	auto m = evaluate(gbn);

	unsigned long long i_max_row = 1;
	unsigned long long i_max_col = 1;
	i_max_col = i_max_col << m->n;
	i_max_row = i_max_row << m->m;

	for(unsigned long long i_row = 0; i_row < i_max_row; i_row++)
		for(unsigned long long i_col = 0; i_col < i_max_col; i_col++)
		{
			if(i_row == i_col)
				REQUIRE(m->get(i_row, i_col) == 1);
			else
				REQUIRE(m->get(i_row, i_col) == 0);
		}
}

TEST_CASE("Evaluating id_1.gbn should yield identity 2 -> 2") {
	auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "id_2.gbn");

	auto m = evaluate(gbn);

	unsigned long long i_max_row = 1;
	unsigned long long i_max_col = 1;
	i_max_col = i_max_col << m->n;
	i_max_row = i_max_row << m->m;

	for(unsigned long long i_row = 0; i_row < i_max_row; i_row++)
		for(unsigned long long i_col = 0; i_col < i_max_col; i_col++)
		{
			if(i_row == i_col)
				REQUIRE(m->get(i_row, i_col) == 1);
			else
				REQUIRE(m->get(i_row, i_col) == 0);
		}
}

TEST_CASE("The evaluation of F-Matrices should result into a diagonal matrix.") {
    auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "f_matrices.gbn");
    auto p_m = evaluate(gbn);
    REQUIRE(p_m->type == DIAGONAL);
}

TEST_CASE("test") {
    auto gbn = read_and_check_gbn(TEST_INSTANCE_FOLDER + "seven_nodes.gbn");

    auto start_time = std::chrono::steady_clock::now();
    auto p_m = evaluate_stepwise(gbn);
    auto end_time = std::chrono::steady_clock::now();
    double diff_milliseconds = std::chrono::duration<double, std::milli>(end_time-start_time).count();
    std::cout << "New: " << diff_milliseconds << std::endl;

    auto start_time2 = std::chrono::steady_clock::now();
    auto p_m2 = evaluate(gbn);
    auto end_time2 = std::chrono::steady_clock::now();
    double diff_milliseconds2 = std::chrono::duration<double, std::milli>(end_time2-start_time2).count();
    std::cout << "Old: " << diff_milliseconds2 << std::endl;


}

TEST_CASE("test1")
{
    std::size_t n_places = 20;
    std::size_t n_transitions = 5;
    std::size_t n_min_tokens = 5;
    std::size_t n_max_tokens = 5;
    std::size_t n_min_pre_places = 1;
    std::size_t n_max_pre_places = 1;
    std::size_t n_min_post_places = 1;
    std::size_t n_max_post_places = 1;

    std::size_t n_simplification_steps = 20;
    std::size_t n_random_transitions_per_simplify = 1;

    std::random_device rd;
    std::mt19937 mt(rd());

    auto cn = randomize_cn(n_places, n_transitions, n_min_tokens, n_max_tokens, n_min_pre_places, n_max_pre_places, n_min_post_places, n_max_post_places, mt);
    auto cn_copy = cn;
    std::vector<std::vector<std::pair<std::size_t, double>>> transitions_w_probabilities;
    std::vector<std::size_t> chosen_transitions;

    auto gbn = build_uniform_independent_obn(n_places);
    auto gbn_copy = gbn;
    auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, 2);

    auto start_time1 = std::chrono::steady_clock::now();
    for(std::size_t i_simplification_step = 0; i_simplification_step < n_simplification_steps; i_simplification_step++)
    {
        for(std::size_t i_rand_transition = 0; i_rand_transition < n_random_transitions_per_simplify; i_rand_transition++)
        {
            auto i_transition = rand_transition_helper.next_p(mt);
            transitions_w_probabilities.push_back(i_transition);

            auto chosen_transition = rand_transition_helper.choose_transition(cn, i_transition);
            chosen_transitions.push_back(chosen_transition);

            fire_with_probability_on_gbn(cn, gbn, i_transition, chosen_transition);
        }
        check_gbn_integrity(gbn);
        simplification(gbn);
        check_gbn_integrity(gbn);
    }
    auto p_m = evaluate(gbn);
    auto end_time1 = std::chrono::steady_clock::now();
    double diff_milliseconds1 = std::chrono::duration<double, std::milli>(end_time1-start_time1).count();
    std::cout << "Old: " << diff_milliseconds1 << std::endl;

    auto start_time = std::chrono::steady_clock::now();
    for(std::size_t i_simplification_step = 0; i_simplification_step < n_simplification_steps; i_simplification_step++) {
        for(std::size_t i_rand_transition = 0; i_rand_transition < n_random_transitions_per_simplify; i_rand_transition++) {
            fire_with_probability_on_gbn(cn_copy, gbn_copy, transitions_w_probabilities.at(i_simplification_step), chosen_transitions.at(i_simplification_step));
        }
        check_gbn_integrity(gbn_copy);
    }
    auto p_m2 = evaluate_stepwise(gbn_copy);
    auto end_time = std::chrono::steady_clock::now();
    double diff_milliseconds = std::chrono::duration<double, std::milli>(end_time-start_time).count();
    std::cout << "New: " << diff_milliseconds << std::endl;

    //Not normalized!!
    //test_matrices_equal(*p_m, *p_m2);
}