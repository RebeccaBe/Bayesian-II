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


/*TEST_CASE("test") {
    std::random_device rd;
    std::mt19937 mt(rd());
    auto cn = randomize_cn(20,10,3,3,1,1,1,1,mt);
    print_cn_details(std::cout, cn);
    std::cout << "Valid transitions: "<< get_n_valid_transitions(cn) << std::endl;
    auto valid_transitions = get_valid_transitions(cn);
    std::cout << "Transitions are: ";
    for(auto tr : valid_transitions)
        std::cout << tr << ", ";

    auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, 3);
    auto list = rand_transition_helper.next_p(mt);
    std::cout << "chosen transitions: ";
    for(auto el : list)
        std::cout << el.first << ", ";
}*/