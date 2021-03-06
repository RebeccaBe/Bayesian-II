#include "../../../libs/catch/catch.hpp"

#include <fstream>

#include "../test_helpers.h"
#include "../../gbn/general/gbn.h"
#include "../../gbn/general/special_cases.h"
#include "../../gbn/evaluation/evaluation.h"
#include "../../gbn/simplification/global_simplification.h"
#include "../../joint_dist/joint_dist.h"
#include "../../joint_dist/special_cases.h"
#include "../../joint_dist/joint_dist_io.h"
#include "../../gbn/matrix/matrix_io.h"
#include "../../gbn/general/gbn_io.h"
#include "../../gbn/general/check.h"
#include "../../gbn/simplification/simplification.h"
#include "../../cn/randomized_generation.h"
#include "../../cn/cn_io.h"
#include "../../cnu/fire_transition.h"
#include "../../cnu/operations_on_gbn.h"
#include "../../cnu/operations_on_joint_dist.h"
#include "../../runtime_tests/helpers/random_transition_helper.h"
#include <random>

TEST_CASE("Uniform dist and uniform gbn should be the same.") {
	auto uniform_joint_dist = build_uniform_joint_dist(3);
	auto uniform_gbn = build_uniform_independent_obn(3);
	auto p_m = evaluate(uniform_gbn);

	test_joint_dist_matrix_equal(uniform_joint_dist, *p_m);
}

TEST_CASE("Joint dist and gbn should be the same 1") {
	JointDist joint_dist;
	for(int i1 = 0; i1 < 2; i1++)
		for(int i2 = 0; i2 < 2; i2++)
			for(int i3 = 0; i3 < 2; i3++)
				if(i1 == 0)
					joint_dist.insert({ { static_cast<bool>(i1), static_cast<bool>(i2), static_cast<bool>(i3) }, 1.0/4 });
				else
					joint_dist.insert({ { static_cast<bool>(i1), static_cast<bool>(i2), static_cast<bool>(i3) }, 0.0 });

	auto gbn = build_independent_obn({{1,0},{0.5,0.5},{0.5,0.5}});
	auto p_m = evaluate(gbn);

	test_joint_dist_matrix_equal(joint_dist, *p_m);
}

TEST_CASE("paper_example.cnu: CNU ops on joint dist and GBN should lead to same dist") 
{
	std::ifstream f_cn(TEST_INSTANCE_FOLDER + "paper_example.cn");
	auto cn = read_cn(f_cn);
	auto cn_copy = cn;
    std::string operation;

    std::ifstream f_joint_dist(TEST_INSTANCE_FOLDER + "paper_example.dist");
	auto dist = read_joint_dist(f_joint_dist);
	std::ifstream f_gbn(TEST_INSTANCE_FOLDER + "paper_example.gbn");
	auto gbn = read_gbn(f_gbn);

	fire_transition_on_gbn(cn, gbn, 3);
	fire_transition_on_joint_dist(cn_copy, dist, 3);

    std::ofstream f("test.dot");
    draw_gbn_graph(f, gbn, "");

    write_gbn(std::cout, gbn);

    simplification(gbn);

	auto p_m = evaluate(gbn);
	print_matrix(std::cout, *p_m);
	print_dist(std::cout, dist);

	test_joint_dist_matrix_equal(dist, *p_m);
}

TEST_CASE("paper_example.cnu: CNU probability ops on joint dist and GBN should lead to same dist")
{
	std::ifstream f_cn(TEST_INSTANCE_FOLDER + "paper_example.cn");
	auto cn = read_cn(f_cn);
	auto cn_copy = cn;

	std::ifstream f_joint_dist(TEST_INSTANCE_FOLDER + "paper_example.dist");
	auto dist = read_joint_dist(f_joint_dist);
	std::ifstream f_gbn(TEST_INSTANCE_FOLDER + "paper_example.gbn");
	auto gbn = read_gbn(f_gbn);

    fire_with_probability_on_gbn(cn, gbn, {{2, 0.75},{1, 0.25}}, 2);
	fire_with_probability_on_joint_dist(cn_copy, dist, {{2, 0.75},{1, 0.25}}, 2);

	simplification(gbn);
	auto p_m = evaluate(gbn);

	test_joint_dist_matrix_equal(dist, *p_m);
}

TEST_CASE("paper_example.cnu: CNU stochastic ops on joint dist and GBN should lead to same dist")
{
    std::ifstream f_cn(TEST_INSTANCE_FOLDER + "paper_example.cn");
    auto cn = read_cn(f_cn);
    auto cn_copy = cn;

    std::ifstream f_joint_dist(TEST_INSTANCE_FOLDER + "paper_example.dist");
    auto dist = read_joint_dist(f_joint_dist);
    std::ifstream f_gbn(TEST_INSTANCE_FOLDER + "paper_example.gbn");
    auto gbn = read_gbn(f_gbn);

    fire_with_probabilityStoch_on_gbn(cn, gbn, {{2, 0.5},{1, 0.2}, {0, 0.3}}, 2);
    fire_with_probabilityStoch_on_joint_dist(cn_copy, dist, {{2, 0.5},{1, 0.2}, {0, 0.3}}, 2);

    simplification(gbn);
    auto p_m = evaluate(gbn);

    test_joint_dist_matrix_equal(dist, *p_m);
}

TEST_CASE("paper_exmple_II.cnu: CNU probabilistic ops on joint dist and GBN should lead to same dist") {

    std::random_device rd;
    std::mt19937 mt(rd());

    std::ifstream f_cn(TEST_INSTANCE_FOLDER + "paper_example_II.cn");
    auto cn = read_cn(f_cn);
    auto cn_copy = cn;

    auto dist = build_uniform_joint_dist(4);

    std::ifstream f_gbn(TEST_INSTANCE_FOLDER + "paper_example_II.gbn");
    auto gbn = read_gbn(f_gbn);

    fire_with_probability_on_gbn(cn_copy, gbn, {{0, 0.25},{1, 0.5}, {2, 0.25}}, 1);
    fire_with_probability_on_joint_dist(cn_copy, dist, {{0, 0.25},{1, 0.5}, {2, 0.25}}, 1);

    auto p_m = evaluate_specific_place(2, gbn);

    test_joint_dist_matrix_equal_marginal_prob(dist, *p_m, 2);
}

TEST_CASE("CNU ops on joint dist and GBN should lead to same dist") 
{
	std::size_t n_places = 22;
	std::size_t n_transitions = 5;
	std::size_t n_min_tokens = 1;
	std::size_t n_max_tokens = 5;
	std::size_t n_min_pre_places = 1;
	std::size_t n_max_pre_places = 2;
	std::size_t n_min_post_places = 1;
	std::size_t n_max_post_places = 2;

	std::size_t n_simplification_steps = 5;
	std::size_t n_random_transitions_per_simplify = 10;

	std::random_device rd;  
	std::mt19937 mt(rd()); 
	std::uniform_int_distribution<std::size_t> rand_transition(0,n_transitions-1);

	auto cn = randomize_cn(n_places, n_transitions, n_min_tokens, n_max_tokens, n_min_pre_places, n_max_pre_places, n_min_post_places, n_max_post_places, mt);
	auto cn_copy = cn;

	auto gbn = build_uniform_independent_obn(n_places);
	auto joint_dist = build_uniform_joint_dist(n_places);

	for(std::size_t i_simplification_step = 0; i_simplification_step < n_simplification_steps; i_simplification_step++)
	{
		for(std::size_t i_rand_transition = 0; i_rand_transition < n_random_transitions_per_simplify; i_rand_transition++)
		{
			auto i_transition = rand_transition(mt);
			fire_transition_on_gbn(cn, gbn, i_transition);
			fire_transition_on_joint_dist(cn_copy, joint_dist, i_transition);
		}
		check_gbn_integrity(gbn);
		simplification(gbn);
		check_gbn_integrity(gbn);
		auto p_m = evaluate(gbn);
		test_joint_dist_matrix_equal(joint_dist, *p_m);
	}
}

TEST_CASE("CNU probability ops on joint dist and GBN should lead to same dist")
{
    std::size_t n_places = 10;
    std::size_t n_transitions = 30;
    std::size_t n_min_tokens = 10;
    std::size_t n_max_tokens = 10;
    std::size_t n_min_pre_places = 1;
    std::size_t n_max_pre_places = 2;
    std::size_t n_min_post_places = 1;
    std::size_t n_max_post_places = 2;

    std::size_t n_simplification_steps = 2;
    std::size_t n_random_transitions_per_simplify = 10;

	std::random_device rd;
	std::mt19937 mt(rd());

	auto cn = randomize_cn(n_places, n_transitions, n_min_tokens, n_max_tokens, n_min_pre_places, n_max_pre_places, n_min_post_places, n_max_post_places, mt);
	auto cn_copy = cn;

	auto gbn = build_uniform_independent_obn(n_places);
	auto joint_dist = build_uniform_joint_dist(n_places);
	auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, 2);
    rand_transition_helper.transition_bubbles = rand_transition_helper.make_transitions_w_probabilities(mt, 0.5);

	for(std::size_t i_simplification_step = 0; i_simplification_step < n_simplification_steps; i_simplification_step++)
	{
		for(std::size_t i_rand_transition = 0; i_rand_transition < n_random_transitions_per_simplify; i_rand_transition++)
		{
			auto i_transition = rand_transition_helper.next_from_bubbles(mt);
			auto chosen_transition = rand_transition_helper.choose_transition(cn, i_transition);
			fire_with_probability_on_gbn(cn, gbn, i_transition, chosen_transition);
			fire_with_probability_on_joint_dist(cn_copy, joint_dist, i_transition, chosen_transition);
		}
		check_gbn_integrity(gbn);
		simplification(gbn);

		check_gbn_integrity(gbn);
		auto p_m = evaluate(gbn);
		test_joint_dist_matrix_equal(joint_dist, *p_m);
	}
}

TEST_CASE("CNU stochastic ops on joint dist and GBN should lead to same dist")
{
    std::size_t n_places = 10;
    std::size_t n_transitions = 30;
    std::size_t n_min_tokens = 5;
    std::size_t n_max_tokens = 10;
    std::size_t n_min_pre_places = 1;
    std::size_t n_max_pre_places = 2;
    std::size_t n_min_post_places = 1;
    std::size_t n_max_post_places = 2;

    std::size_t n_simplification_steps = 2;
    std::size_t n_random_transitions_per_simplify = 10;

    std::random_device rd;
    std::mt19937 mt(rd());

    auto cn = randomize_cn(n_places, n_transitions, n_min_tokens, n_max_tokens, n_min_pre_places, n_max_pre_places, n_min_post_places, n_max_post_places, mt);
    auto cn_copy = cn;

    auto gbn = build_uniform_independent_obn(n_places);
    auto joint_dist = build_uniform_joint_dist(n_places);
    auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, 2);
    rand_transition_helper.transition_bubbles = rand_transition_helper.make_transitions_w_probabilities(mt, 0.5);

    for(std::size_t i_simplification_step = 0; i_simplification_step < n_simplification_steps; i_simplification_step++)
    {
        for(std::size_t i_rand_transition = 0; i_rand_transition < n_random_transitions_per_simplify; i_rand_transition++)
        {
            auto i_transition = rand_transition_helper.next_from_bubbles(mt);
            auto chosen_transition = rand_transition_helper.choose_transition(cn, i_transition);

            fire_with_probabilityStoch_on_gbn(cn, gbn, i_transition, chosen_transition);
            fire_with_probabilityStoch_on_joint_dist(cn_copy, joint_dist, i_transition, chosen_transition);
        }
        check_gbn_integrity(gbn);
        simplification(gbn);
        check_gbn_integrity(gbn);
        auto p_m = evaluate(gbn);

        test_joint_dist_matrix_equal(joint_dist, *p_m);
    }
}

TEST_CASE("Normalizing matrices by columns should work") {
    auto p_matrix = read_matrix({ "dynamic 1 1 [0.2,0;0.2,0]" });

    REQUIRE(p_matrix->get(BitVec("0"), BitVec("0")) == 0.2);
    normalize_matrix_cols(p_matrix);
    REQUIRE(p_matrix->get(BitVec("0"), BitVec("0")) == 0.5);
}

TEST_CASE("CNU ops on joint distribution and GBN should lead to same marginal dist")
{
    std::size_t n_places = 10;
    std::size_t n_transitions = 30;
    std::size_t n_min_tokens = 0;
    std::size_t n_max_tokens = 10;
    std::size_t n_min_pre_places = 1;
    std::size_t n_max_pre_places = 2;
    std::size_t n_min_post_places = 1;
    std::size_t n_max_post_places = 2;

    std::size_t n_simplification_steps = 10;
    std::size_t n_random_transitions_per_simplify = 10;

    std::random_device rd;
    std::mt19937 mt(rd());

    auto cn = randomize_cn(n_places, n_transitions, n_min_tokens, n_max_tokens, n_min_pre_places, n_max_pre_places, n_min_post_places, n_max_post_places, mt);
    auto cn_copy = cn;

    auto gbn = build_uniform_independent_obn(n_places);
    auto joint_dist = build_uniform_joint_dist(n_places);
    auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, 2);
    rand_transition_helper.transition_bubbles = rand_transition_helper.make_transitions_w_probabilities(mt, 0.5);

    for(std::size_t i_simplification_step = 0; i_simplification_step < n_simplification_steps; i_simplification_step++)
    {
        for(std::size_t i_rand_transition = 0; i_rand_transition < n_random_transitions_per_simplify; i_rand_transition++)
        {
            auto i_transition = rand_transition_helper.next_from_bubbles(mt);
            auto chosen_transition = rand_transition_helper.choose_transition(cn, i_transition);
            fire_with_probability_on_gbn(cn, gbn, i_transition, chosen_transition);
            fire_with_probability_on_joint_dist(cn_copy, joint_dist, i_transition, chosen_transition);
        }
        check_gbn_integrity(gbn);

        check_gbn_integrity(gbn);
        auto p_m = evaluate_specific_place(0, gbn);
        auto p_dist = calculate_marginals(0, joint_dist);

        test_joint_dist_matrix_equal_marginal_prob(p_dist, *p_m, 0);
    }
}