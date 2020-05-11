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
#include "../gbn/evaluation/evaluation.h"

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
		("evaluation_type", "Types: DEGREE or FILLIN", cxxopts::value<std::string >()->default_value("degree"))
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

	csv_file << "n_places;milliseconds_degree;milliseconds_fillin" << std::endl;

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

			auto gbn = build_uniform_independent_obn(n_places);
			auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, cn_params.N_MAX_TRANSITIONS_PER_OP);
            rand_transition_helper.transition_bubbles = rand_transition_helper.make_transitions_w_probabilities(mt,0.3);

            auto gbn_copy = gbn;
            auto cn_copy = cn;
            std::vector<std::size_t> chosen_transitions;
            std::vector<std::vector<std::pair<std::size_t, double>>> i_transitions;

            std::string operation;

			auto start_time_degree = std::chrono::steady_clock::now();
			for(std::size_t i_fire = 0; i_fire < cn_params.N_TRANSITIONS_PER_RUN; i_fire++) {
				if(is_detailed)
					std::cout << "i_fire: " << i_fire << std::endl;

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
			evaluate_specific_place(0, gbn, DEGREE);
            auto end_time_degree = std::chrono::steady_clock::now();


            auto start_time_fillIn = std::chrono::steady_clock::now();
            for(std::size_t i_fire = 0; i_fire < cn_params.N_TRANSITIONS_PER_RUN; i_fire++) {
                if(is_detailed)
                    std::cout << "i_fire: " << i_fire << std::endl;

                auto callback = (is_detailed) ? [&operation](std::string high_level, std::string low_level) { std::cout << high_level << " " << low_level << std::endl; } : std::function<void(std::string,std::string)>();
                fire_with_probability_on_gbn(cn_copy, gbn_copy, i_transitions[i_fire], chosen_transitions[i_fire], callback);

                if(is_detailed) {
                    std::ofstream f("run.dot");
                    draw_gbn_graph(f, gbn, std::to_string(i_fire));
                }
            }
            evaluate_specific_place(0, gbn, FILLIN);
            auto end_time_fillIn = std::chrono::steady_clock::now();


			double diff_milliseconds_degree = std::chrono::duration<double, std::milli>(end_time_degree-start_time_degree).count();
            double diff_milliseconds_fillIn = std::chrono::duration<double, std::milli>(end_time_fillIn-start_time_fillIn).count();

			csv_file << n_places << ";" << diff_milliseconds_degree << ";" << diff_milliseconds_fillIn << std::endl;
		}
	}

	return 0;
}