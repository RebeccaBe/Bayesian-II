cmake_minimum_required(VERSION 3.1.0) 

project(Petri_Nets_With_Uncertainty)

# set c++ standard to c++11
set (CMAKE_CXX_STANDARD 17)

# warning level
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -Werror=return-type -Wuninitialized -Wmaybe-uninitialized")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# look for boost library (required)
find_package(Boost REQUIRED regex system filesystem)
include_directories(${Boost_INCLUDE_DIRS})

set(GENERAL_FILES
	src/bitvec.cpp
)

set(CN_FILES
	src/cn/cn_io.cpp
	src/cn/cn_operations.cpp
	src/cn/randomized_generation.cpp
	src/cn/randomized_generation_fc.cpp
        src/cn/cn_stats.cpp)

set(JOINT_DIST_FILES
	src/joint_dist/joint_dist_check.cpp
	src/joint_dist/joint_dist_io.cpp
	src/joint_dist/special_cases.cpp
)

set(CNU_FILES
	src/cnu/operations_on_gbn.cpp
	src/cnu/operations_on_joint_dist.cpp
	src/cnu/fire_transition.cpp
)

set(GBN_FILES
	src/gbn/matrix/matrix.cpp
	src/gbn/matrix/matrix_io.cpp
	src/gbn/matrix/matrix_decomposition.cpp
	src/gbn/matrix/randomized_generation.cpp

	src/gbn/general/gbn.cpp
	src/gbn/general/subgbn.cpp
	src/gbn/general/gbn_io.cpp
	src/gbn/general/check.cpp
	src/gbn/general/path_closing.cpp
	src/gbn/general/randomized_generation.cpp
	src/gbn/general/special_cases.cpp

	src/gbn/evaluation/evaluation.cpp
	src/gbn/evaluation/probability_bookkeeper.cpp
	src/gbn/evaluation/wire_structure.cpp

	src/gbn/modification/vertex_add_remove.cpp
	src/gbn/modification/merging.cpp
	src/gbn/modification/splitting.cpp

	src/gbn/simplification/simplification.cpp
	src/gbn/simplification/local_simplification.cpp
	src/gbn/simplification/global_simplification.cpp
        )

set(RUNTIME_TEST_FILES
	src/runtime_tests/helpers/cn_parameters.cpp
	src/runtime_tests/helpers/random_transition_helper.cpp
)

set(TEST_FILES
	src/tests/cn/cn_io_tests.cpp

	src/tests/cnu/operations_on_joint_dist_tests.cpp
	src/tests/cnu/gbn_and_joint_dist_tests.cpp

	src/tests/gbn/matrix/matrix_io_tests.cpp
	src/tests/gbn/matrix/matrix_decomposition_tests.cpp

	src/tests/gbn/general/gbn_io_tests.cpp

	src/tests/gbn/evaluation/evaluation_tests.cpp
	src/tests/gbn/evaluation/wire_structure_tests.cpp
	src/tests/gbn/evaluation/probability_bookkeeper_tests.cpp

	src/tests/gbn/modification/merging_tests.cpp
	src/tests/gbn/modification/splitting_tests.cpp

        src/tests/gbn/simplification/local_simplification_tests.cpp
	src/tests/gbn/simplification/global_simplification_tests.cpp
	src/tests/gbn/simplification/automated_simplification_tests.cpp
	src/tests/gbn/simplification/automated_local_simplification_tests.cpp
	src/tests/gbn/simplification/automated_global_simplification_tests.cpp

	src/tests/test_helpers.cpp
)

set(TEST_FILES_NEW
	src/tests/cn/cn_fc_tests.cpp

        src/tests/gbn/simplification/local_simplification_tests.cpp

	src/tests/gbn/matrix/matrix_io_tests.cpp

	src/tests/test_helpers.cpp
)

add_library(core_lib STATIC ${GENERAL_FILES} ${CN_FILES} ${JOINT_DIST_FILES} ${CNU_FILES} ${GBN_FILES} ${RUNTIME_TEST_FILES})

add_executable(joint_dist_runner src/main_joint_dist_runner.cpp)
target_link_libraries(joint_dist_runner core_lib)
target_link_libraries(joint_dist_runner Boost::regex Boost::filesystem)

add_executable(gbn_runner src/main_gbn_runner.cpp)
target_link_libraries(gbn_runner core_lib)
target_link_libraries(gbn_runner Boost::regex Boost::filesystem)

add_executable(draw_gbn src/main_draw_gbn.cpp)
target_link_libraries(draw_gbn core_lib)
target_link_libraries(draw_gbn Boost::regex Boost::filesystem)

add_executable(draw_cn src/main_draw_cn.cpp)
target_link_libraries(draw_cn core_lib)
target_link_libraries(draw_cn Boost::regex Boost::filesystem)

add_executable(evaluate_gbn src/main_evaluate_gbn.cpp)
target_link_libraries(evaluate_gbn core_lib)
target_link_libraries(evaluate_gbn Boost::regex Boost::filesystem)

add_executable(merging src/main_merging.cpp)
target_link_libraries(merging core_lib)
target_link_libraries(merging Boost::regex Boost::filesystem)

add_executable(splitting src/main_splitting.cpp)
target_link_libraries(splitting core_lib)
target_link_libraries(splitting Boost::regex Boost::filesystem)

add_executable(simplification src/main_simplification.cpp)
target_link_libraries(simplification core_lib)
target_link_libraries(simplification Boost::regex Boost::filesystem)

add_executable(path_closing src/main_path_closing.cpp)
target_link_libraries(path_closing core_lib)
target_link_libraries(path_closing Boost::regex Boost::filesystem)

add_executable(random_cn src/main_random_cn.cpp)
target_link_libraries(random_cn core_lib)
target_link_libraries(random_cn Boost::regex Boost::filesystem)

add_executable(random_fc_cn src/main_random_fc_cn.cpp)
target_link_libraries(random_fc_cn core_lib)
target_link_libraries(random_fc_cn Boost::regex Boost::filesystem)

add_executable(random_gbn src/main_random_gbn.cpp)
target_link_libraries(random_gbn core_lib)
target_link_libraries(random_gbn Boost::regex Boost::filesystem)

add_executable(sub_gbn src/main_sub_gbn.cpp)
target_link_libraries(sub_gbn core_lib)
target_link_libraries(sub_gbn Boost::regex Boost::filesystem)

add_executable(read_write_gbn src/main_read_write_gbn.cpp)
target_link_libraries(read_write_gbn core_lib)
target_link_libraries(read_write_gbn Boost::regex Boost::filesystem)


add_executable(rt_joint_dist_uniform src/runtime_tests/main_joint_dist_uniform_random_cnu.cpp)
target_link_libraries(rt_joint_dist_uniform core_lib)
target_link_libraries(rt_joint_dist_uniform Boost::regex Boost::filesystem)

add_executable(rt_gbn_uniform src/runtime_tests/main_gbn_uniform_random_cnu.cpp)
target_link_libraries(rt_gbn_uniform core_lib)
target_link_libraries(rt_gbn_uniform Boost::regex Boost::filesystem)

add_executable(rt_demo_simplification src/runtime_tests/main_demo_simplifications.cpp)
target_link_libraries(rt_demo_simplification core_lib)
target_link_libraries(rt_demo_simplification Boost::regex Boost::filesystem)

add_executable(rt_gbn_uniform_p src/runtime_tests/main_gbn_uniform_random_cnu_probability.cpp)
target_link_libraries(rt_gbn_uniform_p core_lib)
target_link_libraries(rt_gbn_uniform_p Boost::regex Boost::filesystem)

add_executable(rt_gbn_uniform_p_marg src/runtime_tests/main_gbn_uniform_random_cnu_probability_marginal.cpp)
target_link_libraries(rt_gbn_uniform_p_marg core_lib)
target_link_libraries(rt_gbn_uniform_p_marg Boost::regex Boost::filesystem)

add_executable(rt_gbn_uniform_p_marg_memory src/runtime_tests/main_gbn_uniform_random_cnu_probability_marginal_matrix_sizes.cpp)
target_link_libraries(rt_gbn_uniform_p_marg_memory core_lib)
target_link_libraries(rt_gbn_uniform_p_marg_memory Boost::regex Boost::filesystem)

add_executable(rt_joint_dist_uniform_p src/runtime_tests/main_joint_dist_uniform_random_cnu_probability.cpp)
target_link_libraries(rt_joint_dist_uniform_p core_lib)
target_link_libraries(rt_joint_dist_uniform_p Boost::regex Boost::filesystem)

add_executable(rt_gbn_uniform_stoch src/runtime_tests/main_gbn_uniform_random_cnu_probabilityStoch.cpp)
target_link_libraries(rt_gbn_uniform_stoch core_lib)
target_link_libraries(rt_gbn_uniform_stoch Boost::regex Boost::filesystem)

add_executable(rt_joint_dist_uniform_stoch src/runtime_tests/main_joint_dist_uniform_random_cnu_probabilityStoch.cpp)
target_link_libraries(rt_joint_dist_uniform_stoch core_lib)
target_link_libraries(rt_joint_dist_uniform_stoch Boost::regex Boost::filesystem)

add_executable(rt_uniform_p src/runtime_tests/main_uniform_probability.cpp)
target_link_libraries(rt_uniform_p core_lib)
target_link_libraries(rt_uniform_p Boost::regex Boost::filesystem)



add_executable(tests src/main_tests.cpp ${TEST_FILES})
target_link_libraries(tests core_lib)
target_link_libraries(tests Boost::regex Boost::filesystem)
target_compile_definitions(tests PRIVATE FOO="${CMAKE_CURRENT_SOURCE_DIR}/src/tests/instances/")

add_executable(tests_new src/main_tests.cpp ${TEST_FILES_NEW})
target_link_libraries(tests_new core_lib)
target_link_libraries(tests_new Boost::regex Boost::filesystem)
target_compile_definitions(tests_new PRIVATE FOO="${CMAKE_CURRENT_SOURCE_DIR}/src/tests/instances/")
