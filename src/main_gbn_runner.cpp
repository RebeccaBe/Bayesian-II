#include <iostream>
#include <fstream>

#include "cn/cn.h"
#include "cn/cn_io.h"
#include "cn/cn_operations.h"
#include "gbn/general/gbn.h"
#include "gbn/general/gbn_io.h"
#include "gbn/general/check.h"
#include "gbn/general/special_cases.h"
#include "gbn/evaluation/evaluation.h"
#include "gbn/simplification/simplification.h"
#include "gbn/matrix/matrix_io.h"
#include "cnu/operations_on_gbn.h"
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <functional>
#include "../libs/cxxopts/include/cxxopts.hpp"
#include "cnu/fire_transition.h"
#include "runtime_tests/helpers/random_transition_helper.h"

// globals (naughty, naughty, ...)
bool is_step_wise_output = false;
bool is_detailed_output = false;
std::string step_wise_path;
std::string detailed_output_folder;
std::size_t n_operations = 0;

void export_gbn(const GBN& gbn, std::string filename)
{
	std::ofstream f(filename);	
	draw_gbn_graph(f, gbn);
}

struct UpdateMetaData {
	std::size_t i_operation;
	std::string high_level_op_string;
	std::string low_level_op_string;
};

void execute_probabilistic_transition(CN &cn, GBN &gbn);

void do_command(std::string command_line, CN& cn, GBN& gbn, std::function<void(UpdateMetaData)> status_callback = std::function<void(UpdateMetaData)>())
{
	boost::trim(command_line);
	if(command_line[0] == 't')
	{
		std::vector<std::string> split_command_line;
		boost::split(split_command_line,command_line,boost::is_any_of(" "));

		auto i_transition = std::stoi(split_command_line[0].substr(1));

		if(status_callback)
			fire_transition_on_gbn(cn, gbn, i_transition, [&status_callback](std::string high_level_op_string, std::string low_level_op_string) { 
				status_callback(UpdateMetaData{ n_operations++, high_level_op_string, low_level_op_string });
			});
		else
			fire_transition_on_gbn(cn, gbn, i_transition);

		return;
	}

	if(command_line.substr(0,4) == "dist")
	{
		auto p_m = evaluate(gbn);
		print_matrix(std::cout, *p_m);
		return;
	}

	if(command_line.substr(0,4) == "draw")
	{
		if(!is_step_wise_output)
			throw std::logic_error("To use 'draw' command please start runner with --step-wise-draw-file");

		std::vector<std::string> tmp_vec;
		boost::split(tmp_vec, command_line, boost::is_any_of(" "));
		std::ofstream f(step_wise_path);
		draw_gbn_graph(f, gbn);
		return;
	}

	if(command_line.substr(0,8) == "simplify" || command_line.substr(0,14) == "simplification")
	{
		simplification(gbn);
		status_callback(UpdateMetaData{ n_operations++, "simplification", "" });
		return;
	}

	if(command_line.substr(0,3) == "cn")
	{
		print_cn_details(std::cout, cn);
		return;
	}

	if(command_line.substr(0,1) == "p")
    {
	    std::cout << "Insert the probability distribution over the transitions as follows: (t0, 2/10) (t1, 3/10) ..." << std::endl
	    << "Note that unmentioned transitions are assigned a probability of 0" << std::endl;

        execute_probabilistic_transition(cn, gbn);
        return;
    }

    if(command_line.substr(0,5) == "stoch")
    {
        std::cout << "Insert the probability distribution over the transitions as follows: (t0, 2/10) (t1, 3/10) ..." << std::endl
                  << "Note that unmentioned transitions are assigned a probability of 0" << std::endl;

        execute_probabilistic_transition(cn, gbn);
        return;
    }

	std::cout << "Command '" << command_line << "' could not be parsed." << std::endl;
}

void execute_probabilistic_transition(CN &cn, GBN &gbn) {
    std::string line;
    std::getline(std::cin, line);
    std::regex regex("^(\\(t([0-9]+), ?([0-9]+)(\\/([0-9])+)?\\) ?)+", std::regex_constants::icase);
    std::smatch matches;
    if(std::regex_match(line, matches, regex)) {
        std::vector<std::pair<std::size_t, double>> transitions_w_probabilities;
        std::vector<std::size_t> used_transitions;

        std::vector<std::string> transitions_strs;
        boost::split(transitions_strs, line, boost::is_any_of(")"));
        transitions_strs.erase(transitions_strs.end());

        for(auto transitions_str : transitions_strs) {
            std::vector<std::string> single_transition_strs;
            boost::trim(transitions_str);

            boost::split(single_transition_strs, transitions_str, boost::is_any_of(","));
            for(auto val_str : single_transition_strs) boost::trim(val_str);
            std::size_t i_transition = std::stoi(single_transition_strs[0].substr(2));

            if(i_transition >= cn.transitions.size()) {
                std::cout << "There were transitions inserted that don't exist. Please insert only transitions between 0 and " << cn.transitions.size()-1 << " ." << std::endl;
                execute_probabilistic_transition(cn, gbn);
                return;
            }

            if(is_in(i_transition, used_transitions)) {
                std::cout << "Some transitions were mentioned multiple times, please try again. To quit, please type 'exit'." << std::endl;
                execute_probabilistic_transition(cn, gbn);
                return;
            }

            transitions_w_probabilities.emplace_back(i_transition, read_double(single_transition_strs[1]));
            used_transitions.push_back(i_transition);
        }

        std::cout << "You can now either specify a transition that gets fired from the available transitions or choose the randomize option by typing 'r'." << std::endl;
        std::cout << "These are the available transitions: ";
        std::vector<std::size_t> valid_transitions;
        for(auto i_transition : transitions_w_probabilities) {
            const auto& transition = cn.transitions[i_transition.first];
            if(check_pre_condition(transition, cn.m)) valid_transitions.push_back(i_transition.first);
        }
        if(valid_transitions.empty()) {
            std::cout << "There are no valid transitions in the given set. Please try again. To quit, please type 'exit'." << std::endl;
            execute_probabilistic_transition(cn, gbn);
            return;
        }
        for(auto i_transition : valid_transitions) std::cout << "t" << i_transition << " ";
        std::cout << std::endl << "If you choose a transition that is not in this list, the randomized option will get executed automatically." << std::endl;
        std::string line2;
        std::getline(std::cin, line2);
        std::size_t chosen_transition = 0;
        if (line2[0] == 't'){
            std::vector<std::string> split_line;
            boost::split(split_line,line2,boost::is_any_of(" "));
            std::size_t transition = std::stoi(split_line[0].substr(1));
            if (is_in(transition, valid_transitions)) {
                chosen_transition = transition;
            } else {
                auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, 0);
                chosen_transition = rand_transition_helper.choose_transition(cn, transitions_w_probabilities);
            }
        } else if(line2[0] == 'r') {
            auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, 0);
            chosen_transition = rand_transition_helper.choose_transition(cn, transitions_w_probabilities);
        } else {
            std::cout << "Command was not understood, please try again. To quit, please type 'exit'." << std::endl;
            execute_probabilistic_transition(cn, gbn);
            return;
        }

            double sum = 0;
        for (auto t : transitions_w_probabilities) sum += t.second;
        if(sum-1 > 0.001 || 1-sum > 0.001) {
            std::cout << "Probabilities do not add up to 1. Instead the sum equals " << std::to_string(sum)
                      << ". Please try again. To quit, please type 'exit'.";
            execute_probabilistic_transition(cn, gbn);
            return;
        }

        std::string operation;
        auto callback = [&operation](std::string high_level, std::string low_level) { std::cout << high_level << " " << low_level << std::endl; };

        fire_with_probability_on_gbn(cn, gbn, transitions_w_probabilities, chosen_transition, callback);
    }
    else if (line.substr(0,4) == "exit") return;
    else {
        std::cout<<"Incorrect syntax, please try again. To quit, please type 'exit'." << std::endl;
        execute_probabilistic_transition(cn, gbn);
    }
}

void execute_stochastic_transition(CN &cn, GBN &gbn) {
    std::string line;
    std::getline(std::cin, line);
    std::regex regex("^(\\(t([0-9]+), ?([0-9]+)(\\/([0-9])+)?\\) ?)+", std::regex_constants::icase);
    std::smatch matches;
    if(std::regex_match(line, matches, regex)) {
        std::vector<std::pair<std::size_t, double>> transitions_w_probabilities;
        std::vector<std::size_t> used_transitions;

        std::vector<std::string> transitions_strs;
        boost::split(transitions_strs, line, boost::is_any_of(")"));
        transitions_strs.erase(transitions_strs.end());

        for(auto transitions_str : transitions_strs) {
            std::vector<std::string> single_transition_strs;
            boost::trim(transitions_str);

            boost::split(single_transition_strs, transitions_str, boost::is_any_of(","));
            for(auto val_str : single_transition_strs) boost::trim(val_str);
            std::size_t i_transition = std::stoi(single_transition_strs[0].substr(2));

            if(i_transition >= cn.transitions.size()) {
                std::cout << "There were transitions inserted that don't exist. Please insert only transitions between 0 and " << cn.transitions.size()-1 << " ." << std::endl;
                execute_probabilistic_transition(cn, gbn);
                return;
            }

            if(is_in(i_transition, used_transitions)) {
                std::cout << "Some transitions were mentioned multiple times, please try again. To quit, please type 'exit'." << std::endl;
                execute_stochastic_transition(cn, gbn);
                return;
            }

            transitions_w_probabilities.emplace_back(i_transition, read_double(single_transition_strs[1]));
            used_transitions.push_back(i_transition);
        }

        std::cout << "You can now either specify a transition that gets fired from the available transitions or choose the randomize option by typing 'r'." << std::endl;
        std::cout << "These are the available transitions: ";
        std::vector<std::size_t> valid_transitions;
        for(auto i_transition : transitions_w_probabilities) {
            const auto& transition = cn.transitions[i_transition.first];
            if(check_pre_condition(transition, cn.m)) valid_transitions.push_back(i_transition.first);
        }
        if(valid_transitions.empty()) {
            std::cout << "There are no valid transitions in the given set. Please try again. To quit, please type 'exit'." << std::endl;
            execute_stochastic_transition(cn, gbn);
            return;
        }
        for(auto i_transition : valid_transitions) std::cout << "t" << i_transition << " ";
        std::cout << std::endl << "If you choose a transition that is not in this list, the randomized option will get executed automatically." << std::endl;
        std::string line2;
        std::getline(std::cin, line2);
        std::size_t chosen_transition = 0;
        if (line2[0] == 't'){
            std::vector<std::string> split_line;
            boost::split(split_line,line2,boost::is_any_of(" "));
            std::size_t transition = std::stoi(split_line[0].substr(1));
            if (is_in(transition, valid_transitions)) {
                chosen_transition = transition;
            } else {
                auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, 0);
                chosen_transition = rand_transition_helper.choose_transition(cn, transitions_w_probabilities);
            }
        } else if(line2[0] == 'r') {
            auto rand_transition_helper = RandomTransitionHelper(cn, RandomTransitionHelper::PROBABILITY, 1, 0);
            chosen_transition = rand_transition_helper.choose_transition(cn, transitions_w_probabilities);
        } else {
            std::cout << "Command was not understood, please try again. To quit, please type 'exit'." << std::endl;
            execute_stochastic_transition(cn, gbn);
            return;
        }

        double sum = 0;
        for (auto t : transitions_w_probabilities) sum += t.second;
        if(sum-1 > 0.001 || 1-sum > 0.001) {
            std::cout << "Probabilities do not add up to 1. Instead the sum equals " << std::to_string(sum)
                      << ". Please try again. To quit, please type 'exit'.";
            execute_stochastic_transition(cn, gbn);
            return;
        }

        std::string operation;
        auto callback = [&operation](std::string high_level, std::string low_level) { std::cout << high_level << " " << low_level << std::endl; };

        fire_with_probabilityStoch_on_gbn(cn, gbn, transitions_w_probabilities, chosen_transition, callback);
    }
    else if (line.substr(0,4) == "exit") return;
    else {
        std::cout<<"Incorrect syntax, please try again. To quit, please type 'exit'." << std::endl;
        execute_stochastic_transition(cn, gbn);
    }
}

CN get_cn(const cxxopts::ParseResult& params)
{
	if(!params.count("cn"))
		throw std::logic_error("Please provide CN file.");

	std::ifstream cn_file(params["cn"].as<std::string>());
	if(!cn_file.is_open())
		throw std::logic_error("CN file could not be opened.");

	return read_cn(cn_file); 
}

GBN get_gbn(const cxxopts::ParseResult& params, const CN& cn)
{
	if(params.count("gbn"))
	{
		std::ifstream gbn_file(params["gbn"].as<std::string>());
		auto gbn = read_gbn(gbn_file);

		if(gbn.n != 0)
			throw std::logic_error("For CNU application, GBN has to have n = 0.");

		if(gbn.m != cn.n)
			throw std::logic_error("For CNU application, GBN has to have m = |cn.places|.");

		return gbn;
	} 

	if(params.count("gbn-uniform-init")) {
		return build_uniform_independent_obn(cn.n);
	}

	throw std::logic_error("Please provide some source of GBN.");
}

int main(int argc, const char** argv)
{
	cxxopts::Options options("random_cn", "Generate random cn instance.");
	options.add_options()
		("help", "Produces this help message.")

		("cn", "Path to .cn file containing petri net model.", cxxopts::value<std::string>())
		("gbn", "Path to .gbn file containing initial gbn.", cxxopts::value<std::string>())
		("gbn-uniform-init", "If this flag is set an initial gbn is used where all vertices are independent and initialiazed with (1/2, 1/2).")
		("step-wise-draw-file", "Path to a .dot file containing the GBN and updated after each 'draw' command.", cxxopts::value<std::string>())
		("detailed-draw-folder", "Path to a folder where after each operation a new .dot file is created.", cxxopts::value<std::string>())
		;
	options.positional_help("<cn> <gbn>").show_positional_help();
	options.parse_positional({"cn", "gbn"});
	auto params = options.parse(argc, argv);

	// show help and exit program if --help is set
	if (params.count("help")) {
		std::cout << options.help() << std::endl;
		return 0;
	}

	auto cn = get_cn(params);
	auto gbn = get_gbn(params,cn);

	is_step_wise_output = params.count("step-wise-draw-file") > 0;
	if(is_step_wise_output) 
		step_wise_path = params["step-wise-draw-file"].as<std::string>();

	is_detailed_output = params.count("detailed-draw-folder") > 0;
	if(is_detailed_output) {
		detailed_output_folder = params["detailed-draw-folder"].as<std::string>();
		if(!detailed_output_folder.empty() && detailed_output_folder.back() != '/')
			detailed_output_folder += '/';
	}

	if(is_detailed_output)
	{
		std::ofstream f(detailed_output_folder + std::to_string(0) + ".dot");
		draw_gbn_graph(f, gbn, "init");
	}

	std::function<void(UpdateMetaData)> callback;

	if(is_detailed_output) {
		callback = [&gbn] (const UpdateMetaData& meta_data) {
			std::cout << meta_data.high_level_op_string << " " << meta_data.low_level_op_string << std::endl;

			std::ofstream f(detailed_output_folder + std::to_string(n_operations) + ".dot");
			draw_gbn_graph(f, gbn, meta_data.high_level_op_string + " " + meta_data.low_level_op_string , true);

			std::ofstream f2(step_wise_path);
			draw_gbn_graph(f2, gbn, meta_data.high_level_op_string + " " + meta_data.low_level_op_string, true);
		};
	} else {
		callback = [&gbn] (const UpdateMetaData& meta_data) {
			std::cout << meta_data.high_level_op_string << " " << meta_data.low_level_op_string << std::endl;
		};
	}

	std::string line;
	while(true) 
	{
		std::getline(std::cin, line);
		do_command(line, cn, gbn, callback);
		std::cout << std::endl;
	}

	return 0;
}
