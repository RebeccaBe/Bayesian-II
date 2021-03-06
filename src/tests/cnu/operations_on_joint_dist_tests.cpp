#include "../../../libs/catch/catch.hpp"

#include <fstream>

#include "../../cn/cn.h"
#include "../../cn/cn_io.h"
#include "../../joint_dist/joint_dist.h"
#include "../../joint_dist/joint_dist_io.h"
#include "../../cnu/operations_on_joint_dist.h"
#include "../test_helpers.h"

TEST_CASE("Set operation should work correctly 1") {
	std::ifstream f1(TEST_INSTANCE_FOLDER + "operations_initial1.dist");
	std::ifstream f2(TEST_INSTANCE_FOLDER + "operations_set1_1.dist");

	auto initial_dist = read_joint_dist(f1); 
	auto expected_dist = read_joint_dist(f2); 

	auto updated_dist = initial_dist;
	set_op({0}, 1, updated_dist);

	compare_joint_dists(updated_dist, expected_dist);
}

TEST_CASE("Set operation should work correctly 2") {
	std::ifstream f1(TEST_INSTANCE_FOLDER + "operations_initial1.dist");
	std::ifstream f2(TEST_INSTANCE_FOLDER + "operations_set1_2.dist");

	auto initial_dist = read_joint_dist(f1); 
	auto expected_dist = read_joint_dist(f2); 

	auto updated_dist = initial_dist;
	set_op({0,2}, 1, updated_dist);

	compare_joint_dists(updated_dist, expected_dist);
}

TEST_CASE("Assert operation should work correctly 1") {
	std::ifstream f1(TEST_INSTANCE_FOLDER + "operations_initial1.dist");
	std::ifstream f2(TEST_INSTANCE_FOLDER + "operations_assert1_1.dist");

	auto initial_dist = read_joint_dist(f1); 
	auto expected_dist = read_joint_dist(f2); 

	auto updated_dist = initial_dist;
	assert_op({0}, 1, updated_dist);

	compare_joint_dists(updated_dist, expected_dist);
}

TEST_CASE("Assert operation should work correctly 2") {
	std::ifstream f1(TEST_INSTANCE_FOLDER + "operations_initial1.dist");
	std::ifstream f2(TEST_INSTANCE_FOLDER + "operations_assert1_2.dist");

	auto initial_dist = read_joint_dist(f1); 
	auto expected_dist = read_joint_dist(f2); 

	auto updated_dist = initial_dist;
	assert_op({0,2}, 0, updated_dist);

	compare_joint_dists(updated_dist, expected_dist);
}

TEST_CASE("Nassert operation should work correctly 1") {
	std::ifstream f1(TEST_INSTANCE_FOLDER + "operations_initial1.dist");
	std::ifstream f2(TEST_INSTANCE_FOLDER + "operations_assert1_1.dist");

	auto initial_dist = read_joint_dist(f1); 
	auto expected_dist = read_joint_dist(f2); 

	auto updated_dist = initial_dist;
	nassert_op({0}, 0, updated_dist);

	compare_joint_dists(updated_dist, expected_dist);
}

TEST_CASE("Nassert operation should work correctly 2") {
	std::ifstream f1(TEST_INSTANCE_FOLDER + "operations_initial1.dist");
	std::ifstream f2(TEST_INSTANCE_FOLDER + "operations_nassert1_1.dist");

	auto initial_dist = read_joint_dist(f1); 
	auto expected_dist = read_joint_dist(f2); 

	auto updated_dist = initial_dist;
	nassert_op({0,2}, 1, updated_dist);

	compare_joint_dists(updated_dist, expected_dist);
}
