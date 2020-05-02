#pragma once

#include "../../cn/cn.h"
#include "../../cn/randomized_generation.h"
#include "../../cn/cn_io.h"
#include "../../cn/cn_operations.h"
#include "../../cn/cn_stats.h"

#include <algorithm>

struct RandomTransitionHelper 
{
	const CN& cn;
	std::vector<std::size_t> i_transitions;
	enum Type {
		UNIFORM,
		FORCED,
		PROBABILITY
	} type;
	double p_success;

	std::size_t max_trans;

	std::size_t n_valid_transitions;

	std::uniform_int_distribution<std::size_t> rand_transition;
	std::uniform_real_distribution<double> uniform_0_1;

    std::vector<std::vector<std::pair<std::size_t, double>>> transition_bubbles;

	RandomTransitionHelper(const CN& cn, Type type, double p_success, std::size_t max_trans = 3);

	std::size_t next(std::mt19937& mt);
	std::vector<std::pair<std::size_t, double>> next_p(std::mt19937& mt);
    std::vector<std::pair<std::size_t, double>> next_from_bubbles(std::mt19937& mt);
    std::size_t choose_transition (CN& cn, std::vector<std::pair<std::size_t, double>> transitions_w_probabilities);

    std::vector<std::vector<std::pair<std::size_t, double>>> make_transitions_w_probabilities(std::mt19937& mt);

};