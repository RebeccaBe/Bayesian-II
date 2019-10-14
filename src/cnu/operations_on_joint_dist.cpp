#include "operations_on_joint_dist.h"
#include <iostream>
#include "../cn/cn.h"
#include "../cn/cn_io.h"
#include "../helpers.hpp"

namespace {
	std::vector<bool> build_target_marking(const std::vector<bool>& m, std::vector<std::size_t> places, bool b)
	{
		auto vec = m;
		for(auto p : places)
			vec[p] = b;

		return vec;
	}

	bool is_equal(const std::vector<bool>& m1, const std::vector<bool>& m2)
	{
		if(m1.size() != m2.size())
			return false;
		for(std::size_t i = 0; i < m1.size(); i++)
			if(m1[i] != m2[i])
				return false;

		return true;
	}
}

void set_op(const std::vector<std::size_t> places, bool b, JointDist& dist)
{
	for(auto& [marking,p] : dist)
	{
		auto target = build_target_marking(marking, places, b);
		if(target != marking)
		{
			dist[target] += p;
			p = 0;
		}
	}
}

void setp_op(const std::size_t place, bool b, JointDist& dist, double probability)
{
	for(auto& [marking,p] : dist)
	{
		auto target_b = build_target_marking(marking, {place}, b);
		if(target_b != marking)
		{
			dist[target_b] += (p*probability);
			p = 0;
		}
		auto target_not_b = build_target_marking(marking, {place}, !b);
		if(target_not_b != marking)
		{
			dist[target_not_b] += (p*(1-probability));
			p = 0;
		}
	}
}

void assert_op(const std::vector<std::size_t> places, bool b, JointDist& dist)
{
	double sum = 0;
	for(auto& [marking,p] : dist)
	{
		auto target = build_target_marking(marking, places, b);
		if(target != marking)
		{
			sum += p;
			p = 0;
		}
	}

	if(sum > 1e-40)
		for(auto& t : dist)
			t.second = t.second/(1-sum);
}

void assert_non_norm_op(const std::vector<std::size_t> places, bool b, JointDist& dist)
{
    for(auto& [marking,p] : dist) {
        auto target = build_target_marking(marking, places, b);
        if(target != marking)
            p = 0;
    }
}

void normalize_op(JointDist& dist)
{
    double sum = 0;
    for(auto& [marking,p] : dist)
        sum += p;

    if(sum > 1e-40)
        for(auto& t : dist)
            t.second = t.second/sum;
}

void nassert_op(const std::vector<std::size_t> places, bool b, JointDist& dist)
{
	double sum = 0;
	for(auto& [marking,p] : dist)
	{
		auto target = build_target_marking(marking, places, b);
		if(target == marking)
		{
			sum += p;
			p = 0;
		}
	}

	if(sum > 1e-40)
		for(auto& t : dist)
			t.second = t.second/(1-sum);
}

void successp_op(const std::vector<std::vector<std::size_t>> pre_places, const std::vector<std::vector<std::size_t>> post_places,
				const std::vector<double> probabilities, JointDist& dist) {

    double sum = 0;
    for(auto p : probabilities) sum += p;
    if(sum-1 > 0.001 || 1-sum > 0.001)
        throw std::logic_error(std::string("Probabilities do not add up to 1. Instead the sum equals " + std::to_string(sum)));


    auto original_dist = dist;
	for (auto& [marking,p] : dist) {dist[marking] = 0;}

	for(std::size_t i = 0; i < probabilities.size(); i++) {
		double probability = probabilities[i];

		auto worker_dist = original_dist;
		assert_non_norm_op(pre_places[i], 1, worker_dist);
		set_op(pre_places[i], 0, worker_dist);
		set_op(post_places[i], 1, worker_dist);

		for (auto& [marking,p] : dist) p += (worker_dist.at(marking) * probability);
	}

	normalize_op(dist);
}
