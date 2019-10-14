#include "../../../libs/catch/catch.hpp"

#include <fstream>
#include <iostream>
#include <chrono>

#include "../../cn/cn.h"
#include "../../cn/cn_io.h"
#include "../../cn/randomized_generation_fc.h"
#include "../../cn/randomized_generation.h"
#include "../../cn/cn_stats.h"
#include "../../helpers.hpp"
#include "../test_helpers.h"
#include "../../runtime_tests/helpers/random_transition_helper.h"

#include "../../gbn/general/gbn.h"
#include "../../gbn/general/gbn_io.h"
#include "../../gbn/evaluation/evaluation.h"
#include "../../gbn/matrix/matrix_io.h"
#include "../../gbn/matrix/matrix.h"
#include "../../gbn/simplification/simplification.h"
#include "../../gbn/general/subgbn.h"
#include "../../gbn/modification/merging.h"

#include "../../joint_dist/joint_dist.h"
#include "../../joint_dist/special_cases.h"
#include "../../joint_dist/joint_dist_io.h"

#include "../../cnu/fire_transition.h"
#include "../../cnu/operations_on_gbn.h"
#include "../../cnu/operations_on_joint_dist.h"
#include "../../gbn/general/special_cases.h"

TEST_CASE("Updating transitions in a free choice net should work.") {

    CN initial_cn;
    Transition t1, t2, t3, t4;
    t1.pre = std::vector<std::size_t>{1, 2, 3};
    t1.post = std::vector<std::size_t>{4};

    t2.pre = std::vector<std::size_t>{7, 9};
    t2.post = std::vector<std::size_t>{4};

    t3.pre = std::vector<std::size_t>{5, 6};
    t3.post = std::vector<std::size_t>{4};

    t4.pre = std::vector<std::size_t>{1, 2, 3};
    t4.post = std::vector<std::size_t>{4};

    initial_cn.transitions.push_back(t1);
    initial_cn.transitions.push_back(t2);
    initial_cn.transitions.push_back(t3);
    initial_cn.transitions.push_back(t4);

    std::vector<std::size_t> places{1, 2, 3, 7, 9};
    std::vector<std::size_t> new_result1{1, 2, 3, 7, 9};
    std::vector<std::size_t> new_result2{7, 9, 1, 2, 3};

    auto updated_cn = initial_cn;
    update_transitions(updated_cn, places);

    REQUIRE(updated_cn.transitions.size() == 4);
    REQUIRE(updated_cn.transitions[0].pre == new_result1);
    REQUIRE(updated_cn.transitions[1].pre == new_result1);
    REQUIRE(updated_cn.transitions[2].pre == initial_cn.transitions[2].pre);
    REQUIRE(updated_cn.transitions[3].pre == new_result1);
    REQUIRE(places == new_result1);
}

/*TEST_CASE("generation free choice net") {
    std::random_device rd;
    std::mt19937 mt(rd());

    auto cn = randomize_fc_cn(20, 30, 5, 15, 2, 3, 2, 3, mt);
    auto cn1 = cn;
    eliminate_loops(cn1);
    std::cout << "valid transitions: " << get_n_valid_transitions(cn)  << std::endl;
    std::cout << "transitions: " << get_n_transitions(cn) << std::endl;
    std::cout << "number of markings: " << get_n_markings(cn) << std::endl;
    auto pre_list = get_distribution_pre(cn);
    auto post_list = get_distribution_post(cn);
    for(auto pre:pre_list)
        std::cout << pre << ", ";
    std::cout << std::endl;
    for(auto post:post_list)
        std::cout << post << ", ";
    std::cout << "\n average of degree of pre: " << get_average_pre(cn) << std::endl;
    std::cout << "average of degree of post: " << get_average_post(cn) << std::endl;

    std::ofstream f("free_choice.txt");
    std::ofstream f1("free_choice_updated.txt");
    export_cn(f, cn);
    export_cn(f1, cn1);
}*/

TEST_CASE("testing time difference") {
    std::random_device rd;
    std::mt19937 mt(rd());

    std::ofstream csv_file("time_testing.csv");

    csv_file << "n_places,milli normal,valid transitions" << std::endl;

    for(std::size_t n_places = 10; n_places <= 20; n_places++) {
        for(std::size_t i_run = 0; i_run < 10; i_run++) {
            std::size_t n_transitions = std::llround(n_places * 5);
            std::size_t n_min_tokens = std::llround(n_places * 0.1);
            std::size_t n_max_tokens = std::llround(n_places * 0.9);

            auto start_time = std::chrono::steady_clock::now();

            auto cn = randomize_cn(n_places, n_transitions, n_min_tokens, n_max_tokens, 1, 3, 1, 3, mt);

            auto end_time = std::chrono::steady_clock::now();
            double diff_milliseconds = std::chrono::duration<double, std::milli>(end_time - start_time).count();

            csv_file << n_places << "," << diff_milliseconds << "," << get_n_valid_transitions(cn) << std::endl;
        }
    }

    csv_file << "n_places,milli free choice,valid transitions" << std::endl;

    for(std::size_t n_places = 10; n_places <= 20; n_places++) {
        for(std::size_t i_run = 0; i_run < 10; i_run++) {
            std::size_t n_transitions = std::llround(n_places * 5);
            std::size_t n_min_tokens = std::llround(n_places * 0.1);
            std::size_t n_max_tokens = std::llround(n_places * 0.9);

            auto start_time = std::chrono::steady_clock::now();

            auto cn = randomize_fc_cn(n_places, n_transitions, n_min_tokens, n_max_tokens, 1, 3, 1, 3, mt);

            auto end_time = std::chrono::steady_clock::now();
            double diff_milliseconds = std::chrono::duration<double, std::milli>(end_time - start_time).count();

            csv_file << n_places << "," << diff_milliseconds << "," << get_n_valid_transitions(cn) << std::endl;
        }
    }
}

TEST_CASE("test") {
    std::ifstream f(TEST_INSTANCE_FOLDER + "uncertainty_test.cnu");
    auto cn = read_cn(f);

    auto gbn = build_uniform_independent_obn(cn.n);

    fire_transition_on_gbn(cn, gbn, 0);
    fire_transition_on_gbn(cn, gbn, 1);
    write_gbn(std::cout, gbn);


    successp_op({{0},{4},{2}}, {{1},{5},{3}}, {0.2,0.7,.1}, gbn);
    write_gbn(std::cout, gbn);
}