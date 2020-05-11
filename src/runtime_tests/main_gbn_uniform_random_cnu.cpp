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
#include "../cn/randomized_generation_fc.h"
#include "../cn/randomized_generation.h"
#include "../cn/cn_io.h"
#include "../cn/cn_operations.h"
#include "../cn/cn_stats.h"
#include "../cnu/fire_transition.h"
#include "helpers/random_transition_helper.h"
#include "helpers/cn_parameters.h"

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
		("rand-transition-type", "", cxxopts::value<std::string>()->default_value("uniform"))
		("n-max-transitions-per-op", "", cxxopts::value<std::size_t>()->default_value("3"))

		("export-name", "", cxxopts::value<std::string>()->default_value("out"))
		("detailed", "")

		("evaluation_type", "Type cannot be changed", cxxopts::value<std::string>()->default_value("default"))
		("free-choice", "1 or 0", cxxopts::value<bool>()->default_value("false"))
		//("comparison-value", "Compared to n = how many transition get set off once at a time?", cxxopts::value<std::size_t>()->default_value("1"))
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
    std::map<std::string, std::size_t> operations;
    std::map<std::string, std::size_t> operations2;

	std::ofstream csv_file(params["export-name"].as<std::string>()+".csv");

	csv_file << "n_places;milliseconds;valid_transitions" << std::endl;
    //csv_file << "n_places;milliseconds;milliseconds_n_at_once;valid_transitions" << std::endl;

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

			auto valid_transitions = get_n_valid_transitions(cn);
			auto gbn = build_uniform_independent_obn(n_places);
			auto rand_transition_helper = RandomTransitionHelper(cn, cn_params.RAND_TRANSITION_TYPE, cn_params.P_SUCCESS);

            std::vector<std::size_t> transition_order;
            auto cn2 = cn;
            auto gbn2 = gbn;
            std::string operation;

			auto start_time = std::chrono::steady_clock::now();

			for(std::size_t i_fire = 0; i_fire < cn_params.N_TRANSITIONS_PER_RUN; i_fire++)
			{
				if(is_detailed)
					std::cout << "i_fire: " << i_fire << std::endl;
				auto i_transition = rand_transition_helper.next(mt);
				//transition_order.push_back(i_transition);
				auto callback = (is_detailed) ? [&operation](std::string high_level, std::string low_level) { std::cout << high_level << " " << low_level << std::endl; } : std::function<void(std::string,std::string)>();
                auto callback_simplification = (is_detailed) ? [&operations](const GBN& gbn,std::string op) { operations[op]++;} : std::function<void(const GBN&,std::string)>();
				fire_transition_on_gbn(cn, gbn, i_transition, callback);

				old_simplification(gbn, callback_simplification);

				if(is_detailed)
				{
					std::ofstream f("run.dot");
					draw_gbn_graph(f, gbn, std::to_string(i_fire));
				}
			}
            auto end_time = std::chrono::steady_clock::now();


           /* auto start_time2 = std::chrono::steady_clock::now();

            for(std::size_t i_fire = 0; i_fire < cn_params.N_TRANSITIONS_PER_RUN; i_fire++)
            {

                if(is_detailed)
                    std::cout << "i_fire: " << i_fire << std::endl;
                auto i_transition = transition_order[i_fire];
                auto callback = (is_detailed) ? [&operation](std::string high_level, std::string low_level) { std::cout << high_level << " " << low_level << std::endl; operation = high_level; } : std::function<void(std::string,std::string)>();
                auto callback_simplification = (is_detailed) ? [&operations2](const GBN& gbn,std::string op) { operations2[op]++;} : std::function<void(const GBN&,std::string)>();
                fire_transition_on_gbn(cn2, gbn2, i_transition, callback);

                if((i_fire != 0 && (i_fire+1) % cn_params.COMPARISON_VALUE == 0) || i_fire == (cn_params.N_TRANSITIONS_PER_RUN-1))
                    simplification(gbn2, callback_simplification);

                if(is_detailed)
                {
                    std::ofstream f("runBlock.dot");
                    draw_gbn_graph(f, gbn2, std::to_string(i_fire));
                }
            }
            auto end_time2 = std::chrono::steady_clock::now();*/

			double diff_milliseconds = std::chrono::duration<double, std::milli>(end_time-start_time).count();
            //double diff_milliseconds2 = std::chrono::duration<double, std::milli>(end_time2-start_time2).count();

			csv_file << n_places << ";" << diff_milliseconds << ";" //<< diff_milliseconds2 << ";"
			    << valid_transitions << std::endl;

		}
	}

	if(is_detailed) {
        for (auto operation : operations) std::cout << operation.first << ": " << operation.second << std::endl;
        std::cout << "block firing: " << std::endl;
        for (auto operation : operations2) std::cout << operation.first << ": " << operation.second << std::endl;
    }

	return 0;
}