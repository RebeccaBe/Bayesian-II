#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <iomanip>

#include "../../libs/cxxopts/include/cxxopts.hpp"

#include "../gbn/general/gbn.h"
#include "../gbn/general/gbn_io.h"
#include "../gbn/general/special_cases.h"
#include "../gbn/simplification/simplification.h"
#include "../gbn/evaluation/evaluation.h"
#include "../joint_dist/joint_dist.h"
#include "../joint_dist/special_cases.h"
#include "../joint_dist/joint_dist_io.h"
#include "../cn/randomized_generation_fc.h"
#include "../cn/randomized_generation.h"
#include "../cn/cn_io.h"
#include "../cn/cn_operations.h"
#include "../cn/cn_stats.h"
#include "../cnu/fire_transition.h"
#include "helpers/random_transition_helper.h"
#include "helpers/cn_parameters.h"

namespace {
    bool bool_vector_reverse_compare(const std::vector<bool>& v1, const std::vector<bool>& v2)
    {
        std::size_t min_size = std::min(v1.size(), v2.size());
        for(long long i = min_size-1; i >= 0; i--)
        {
            if(v1[i] != v2[i])
                return v1[i] < v2[i];
        }

        return v1.size() < v2.size();
    }
}

bool test_joint_dist_matrix_equal_marginal_prob(const JointDist& joint_dist, const Matrix& m1, std::size_t place)
{
    bool equal = true;

    if(m1.n != 0)
        throw std::logic_error("Cannot compare joint dist with matrix with n != 0");

    std::vector<std::pair<std::vector<bool>, double>> joint_dist_ordered;
    for(auto t : joint_dist)
        joint_dist_ordered.push_back(t);

    std::sort(joint_dist_ordered.begin(), joint_dist_ordered.end(), [](const auto& p1, const auto& p2) { return bool_vector_reverse_compare(p1.first, p2.first); });

    unsigned long long to_max_dist = joint_dist_ordered.size();
    std::size_t matrix_max = (log10(to_max_dist))/(log10(2));
    MatrixPtr matrix = std::make_shared<DynamicMatrix>(0,matrix_max);

    for(Index to = 0; to < to_max_dist; to++) {
        matrix->set(to, 0, joint_dist_ordered[to].second);
    }

    Matrix& m2 = *matrix;

    if(m1.n > 0 || m2.n > 0)
        throw std::logic_error("Both matrices cannot have more than one column (n=0).");

    MatrixPtr sum_matrix = std::make_shared<DynamicMatrix>(0, 1);

    unsigned long long to_max_sum = 1;
    to_max_sum = to_max_sum << m2.m;

    for(Index to = 0; to < to_max_sum; to++) {
        BitVec assignment = to;
        if(assignment[place] == 0)
            sum_matrix->add(0,0,m2.get(to, 0));
        else if (assignment[place] == 1)
            sum_matrix->add(1,0,m2.get(to, 0));
    }

    Matrix& sum_m = *sum_matrix;

    unsigned long long to_max = 1;
    unsigned long long from_max = 1;
    to_max = to_max << m1.m;
    from_max = from_max << m1.n;
    for(Index to = 0; to < to_max; to++)
        for(Index from = 0; from < from_max; from++)
            if(!(m1.get(to, from) - (sum_m.get(to, from)) < 0.0001) || !((sum_m.get(to, from)) - m1.get(to, from) < 0.0001))
                equal = false;

    return equal;
}

std::size_t max_matrix_dimension(GBN gbn) {
    auto g = gbn.graph;
    std::size_t max = 0;
    for(auto v : all_vertices(gbn)) {
        if(type(v, g) != OUTPUT) {
            auto m = matrix(v, gbn.graph);
            switch(m->type) {
                case DIAGONAL:
                case F:
                    if (m->n > max)
                        max = m->n;
                    break;
                default:
                    if ((m->n)+(m->m) > max)
                        max = (m->n)+(m->m);
            }
        }
    }
    return max;
}

int main(int argc, const char** argv)
{
    cxxopts::Options options("random_cn", "Generate random cn instance.");
    options.add_options()
            ("help", "Produces this help message.")

            ("n-min-places", "", cxxopts::value<std::size_t>()->default_value("10"))
            ("n-max-places", "", cxxopts::value<std::size_t>()->default_value("20"))

            ("n-runs", "", cxxopts::value<std::size_t>()->default_value("10"))
            ("n-transitions-per-run", "", cxxopts::value<std::size_t>()->default_value("100"))

            ("percent-min-tokens", "", cxxopts::value<double>()->default_value("0.1"))
            ("percent-max-tokens", "", cxxopts::value<double>()->default_value("0.9"))

            ("percent-transitions", "", cxxopts::value<std::size_t>()->default_value("5"))

            ("n-min-pre-places", "", cxxopts::value<std::size_t>()->default_value("1"))
            ("n-max-pre-places", "", cxxopts::value<std::size_t>()->default_value("3"))
            ("n-min-post-places", "", cxxopts::value<std::size_t>()->default_value("1"))
            ("n-max-post-places", "", cxxopts::value<std::size_t>()->default_value("3"))

            ("p-success", "", cxxopts::value<double>()->default_value("0.5"))
            ("rand-transition-type", "Type cannot be changed.", cxxopts::value<std::string>()->default_value("probability"))
            ("n-max-transitions-per-op", "", cxxopts::value<std::size_t>()->default_value("3"))

            ("export-name", "", cxxopts::value<std::string>()->default_value("out"))
            ("detailed", "")

            ("free-choice", "1 or 0", cxxopts::value<bool>()->default_value("false"))
            ("evaluation_type", "Types: degree or fillIn", cxxopts::value<std::string>()->default_value("degree"))
            ;

    auto params = options.parse(argc, argv);

    if (params.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    auto cn_params = read_cn_params(params);
    std::ofstream param_file(params["export-name"].as<std::string>()+".params");
    write_cn_params(param_file,cn_params);

    bool is_detailed = params.count("detailed") > 0;
    bool is_free_choice = cn_params.FREE_CHOICE;

    std::random_device rd;
    std::mt19937 mt(rd());

    std::ofstream csv_file(params["export-name"].as<std::string>()+".csv");

    csv_file << "n_places;milliseconds_gbn;milliseconds_dist;max_dimension" << std::endl;

    for(std::size_t n_places = cn_params.N_MIN_PLACES; n_places <= cn_params.N_MAX_PLACES; n_places++)
    {
        if(is_detailed)
            std::cout << "n_places: " << n_places << std::endl;
        for(std::size_t i_run = 0; i_run < cn_params.N_RUNS; i_run++)
        {
            if(is_detailed)
                std::cout << "i_run: " << i_run << std::endl;
            std::size_t n_min_tokens = std::llround(n_places * cn_params.PERCENT_MIN_TOKENS);
            std::size_t n_max_tokens = std::llround(n_places * cn_params.PERCENT_MAX_TOKENS);
            std::size_t n_transitions = std::llround(n_places * cn_params.PERCENT_TRANSITIONS);

            // generate random cn
            CN cn;
            if(is_free_choice) {
                cn = randomize_fc_cn(n_places, n_transitions, n_min_tokens, n_max_tokens, cn_params.N_MIN_PRE_PLACES, cn_params.N_MAX_PRE_PLACES, cn_params.N_MIN_POST_PLACES, cn_params.N_MAX_POST_PLACES, mt);
            } else {
                cn = randomize_cn(n_places, n_transitions, n_min_tokens, n_max_tokens, cn_params.N_MIN_PRE_PLACES, cn_params.N_MAX_PRE_PLACES, cn_params.N_MIN_POST_PLACES, cn_params.N_MAX_POST_PLACES, mt);
            }
            auto cn_copy = cn;
            auto cn_copy_copy = cn;

            std::vector<std::size_t> chosen_transitions;
            std::vector<std::vector<std::pair<std::size_t, double>>> i_transitions;

            auto gbn = build_uniform_independent_obn(n_places);
            auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, cn_params.N_MAX_TRANSITIONS_PER_OP);
            rand_transition_helper.transition_bubbles = rand_transition_helper.make_transitions_w_probabilities(mt,1);

            std::string operation;

            auto start_time_gbn = std::chrono::steady_clock::now();
            for(std::size_t i_fire = 0; i_fire < cn_params.N_TRANSITIONS_PER_RUN; i_fire++) {
                if(is_detailed)
                    std::cout << "gbn i_fire: " << i_fire << std::endl;

                auto i_transition = rand_transition_helper.next_from_bubbles(mt);
                auto chosen_transition = rand_transition_helper.choose_transition(cn, i_transition);

                i_transitions.push_back(i_transition);
                chosen_transitions.push_back(chosen_transition);

                auto callback = (is_detailed) ? [&operation](std::string high_level, std::string low_level) { std::cout << high_level << " " << low_level << std::endl; } : std::function<void(std::string,std::string)>();
                fire_with_probability_on_gbn(cn, gbn, i_transition, chosen_transition, callback);

                if(is_detailed) {
                    std::ofstream f("run.dot");
                    draw_gbn_graph(f, gbn, std::to_string(i_fire));
                }
            }
            auto p_m = evaluate_specific_place(0, gbn, cn_params.EVALUATION_TYPE);
            auto end_time_gbn = std::chrono::steady_clock::now();

            auto joint_dist = build_uniform_joint_dist(n_places);
            auto start_time_dist = std::chrono::steady_clock::now();
            for(std::size_t i_fire = 0; i_fire < cn_params.N_TRANSITIONS_PER_RUN; i_fire++) {
                if(is_detailed)
                    std::cout << "dist i_fire: " << i_fire << std::endl;

                auto callback = (is_detailed) ? [&operation](std::string high_level, std::string low_level) { std::cout << high_level << " " << low_level << std::endl; } : std::function<void(std::string,std::string)>();
                fire_with_probability_on_joint_dist(cn_copy_copy, joint_dist, i_transitions[i_fire], chosen_transitions[i_fire], callback);
            }
            calculate_marginals(0, joint_dist);
            auto end_time_dist = std::chrono::steady_clock::now();


            double diff_milliseconds_gbn = std::chrono::duration<double, std::milli>(end_time_gbn-start_time_gbn).count();
            double diff_milliseconds_dist = std::chrono::duration<double, std::milli>(end_time_dist-start_time_dist).count();


            csv_file << n_places << ";" << diff_milliseconds_gbn << ";" << diff_milliseconds_dist << ";" << std::to_string(max_matrix_dimension(gbn)) << std::endl;

            if(is_detailed) {
                if(!test_joint_dist_matrix_equal_marginal_prob(joint_dist, *p_m, 0))
                    std::cout << "Matrices are not equal! (degree)" << std::endl;
            }
        }
    }

    return 0;
}