#include "fire_transition.h"

#include "../cn/cn_operations.h"
#include "../cnu/operations_on_gbn.h"
#include "../cnu/operations_on_joint_dist.h"
#include "../gbn/general/check.h"
#include <iostream>
#include <random>

void fire_transition_on_gbn(CN& cn, GBN& gbn, std::size_t i_transition, std::function<void(std::string,std::string)> status_callback)
{
	const auto& transition = cn.transitions[i_transition];
	const auto& pre_places = transition.pre;
	const auto& post_places = transition.post;

	// fail_pre case
	if(!check_pre_condition(transition,cn.m))
	{
		nassert_op(transition.pre, 1, gbn);
		check_gbn_integrity(gbn);

		if(status_callback)
			status_callback(std::string("fail-pre_{t")+std::to_string(i_transition)+"}", std::string("nassert"));
		return;
	}

	assert_op(pre_places,1, gbn);
	check_gbn_integrity(gbn);
	if(status_callback)
		status_callback(std::string("success_{t")+std::to_string(i_transition)+"}", std::string("assert_pre"));

	set_op(pre_places,0,gbn);
	check_gbn_integrity(gbn);
	if(status_callback)
		status_callback(std::string("success_{t")+std::to_string(i_transition)+"}", std::string("set_pre"));

	set_op(post_places,1,gbn);
	check_gbn_integrity(gbn);
	if(status_callback)
		status_callback(std::string("success_{t")+std::to_string(i_transition)+"}", std::string("set_post"));

	// update marking in CNU
	for(const auto& p : transition.pre)
		cn.m.at(p) = 0;
	for(const auto& p : transition.post)
		cn.m.at(p) = 1;
}

void fire_transition_on_joint_dist(CN& cn, JointDist& dist, std::size_t i_transition, std::function<void(std::string,std::string)> status_callback)
{
	const auto& transition = cn.transitions[i_transition];
	const auto& pre_places = transition.pre;
	const auto& post_places = transition.post;

	// fail_pre case
	if(!check_pre_condition(transition,cn.m))
	{
		nassert_op(transition.pre, 1, dist);

		if(status_callback)
			status_callback(std::string("fail-pre_{t")+std::to_string(i_transition)+"}", std::string("nassert"));
		return;
	}

	assert_op(pre_places,1, dist);
	if(status_callback)
		status_callback(std::string("success_{t")+std::to_string(i_transition)+"}", std::string("assert_pre"));

	set_op(pre_places,0,dist);
	if(status_callback)
		status_callback(std::string("success_{t")+std::to_string(i_transition)+"}", std::string("set_pre"));

	set_op(post_places,1,dist);
	if(status_callback)
		status_callback(std::string("success_{t")+std::to_string(i_transition)+"}", std::string("set_post"));

	// update marking in CNU
	for(const auto& p : transition.pre)
		cn.m[p] = 0;
	for(const auto& p : transition.post)
		cn.m[p] = 1;
}

void fire_with_probability_on_gbn (CN& cn, GBN& gbn, std::vector<std::pair<std::size_t, double>> transitions, std::size_t chosen_transition, std::function<void(std::string,std::string)> status_callback) {
    if(transitions.empty()) {
        std::cout << "probability firing: No transitions have been chosen." << std::endl;
        return;
    }

    std::vector<std::pair<std::size_t, double>> valid_transitions;
    for(auto i_transition : transitions) {
        const auto& transition = cn.transitions[i_transition.first];
        if(check_pre_condition(transition, cn.m)) valid_transitions.push_back(i_transition);
    }

    if(valid_transitions.empty()) {

        if(transitions.size() == 1) {
            const auto& transition = cn.transitions[transitions[0].first];
            const auto& pre_places = transition.pre;

            nassert_op(pre_places, true, gbn);
            check_gbn_integrity(gbn);
            if(status_callback)
                status_callback(std::string("fail_{t")+std::to_string(transitions[0].first)+"}", std::string("nassert_pre"));

            return;
        }

        std::vector<double> probabilities;
        std::vector<std::vector<std::size_t>> all_pre_places;
        std::string transitions_string = "";

        for(auto [i_transition, probability] : transitions) {
            const auto transition = cn.transitions[i_transition];

            probabilities.push_back(probability);
            all_pre_places.push_back(transition.pre);
            if(status_callback) {
                std::stringstream ss;
                ss << "t" << i_transition << ": " << probability << ", ";
                transitions_string.append(ss.str());
            }
        }

        auto& pre_places = all_pre_places;
        failp_op(pre_places, probabilities, gbn);
        check_gbn_integrity(gbn);
        if(status_callback) status_callback(std::string("fail_p{")+ transitions_string +"}", std::string("fail_p"));

        return;
    }

    if(transitions.size() == 1) {
		const auto& transition = cn.transitions[transitions[0].first];
		const auto& pre_places = transition.pre;
		const auto& post_places = transition.post;

		assert_op(pre_places,1, gbn);
		check_gbn_integrity(gbn);
		if(status_callback)
			status_callback(std::string("success_{t")+std::to_string(transitions[0].first)+"}", std::string("assert_pre"));

		set_op(pre_places,0,gbn);
		check_gbn_integrity(gbn);
		if(status_callback)
			status_callback(std::string("success_{t")+std::to_string(transitions[0].first)+"}", std::string("set_pre"));

		set_op(post_places,1,gbn);
		check_gbn_integrity(gbn);
		if(status_callback)
			status_callback(std::string("success_{t")+std::to_string(transitions[0].first)+"}", std::string("set_post"));

		if(valid_transitions.empty())
			return;

		// update marking in CNU
		for(const auto& p : transition.pre)
			cn.m.at(p) = 0;
		for(const auto& p : transition.post)
			cn.m.at(p) = 1;
    	return;
	}

	std::vector<double> probabilities;
	std::vector<std::vector<std::size_t>> all_pre_places;
	std::vector<std::vector<std::size_t>> all_post_places;
	std::string transitions_string = "";

	for(auto [i_transition, probability] : transitions) {
		const auto transition = cn.transitions[i_transition];

		probabilities.push_back(probability);
		all_pre_places.push_back(transition.pre);
		all_post_places.push_back(transition.post);
		if(status_callback) {
			std::stringstream ss;
			ss << "t" << i_transition << ": " << probability << ", ";
			transitions_string.append(ss.str());
		}
	}

	auto& pre_places = all_pre_places;
	auto& post_places = all_post_places;
	successp_op(pre_places, post_places, probabilities, gbn);
	check_gbn_integrity(gbn);
	if(status_callback)
		status_callback(std::string("successp_{")+ transitions_string +"}", std::string("successp"));

    if(valid_transitions.empty())
        return;

	// update marking in CNU
	if(!valid_transitions.empty()){
        const auto& transition = cn.transitions[chosen_transition];
        for(const auto& p : transition.pre)
            cn.m.at(p) = 0;
        for(const auto& p : transition.post)
            cn.m.at(p) = 1;
	}
}

void fire_with_probability_on_joint_dist(CN& cn, JointDist& dist, std::vector<std::pair<std::size_t, double>> transitions, std::size_t chosen_transition, std::function<void(std::string,std::string)> status_callback) {
    if(transitions.empty()) {
        std::cout << "probability firing: No transitions have been chosen." << std::endl;
        return;
    }

    std::vector<std::pair<std::size_t, double>> valid_transitions;
    for(auto i_transition : transitions) {
        const auto& transition = cn.transitions[i_transition.first];
        if(check_pre_condition(transition, cn.m)) valid_transitions.push_back(i_transition);
    }

    if(valid_transitions.empty()) {

        if(transitions.size() == 1) {
            const auto& transition = cn.transitions[transitions[0].first];
            const auto& pre_places = transition.pre;

            nassert_op(pre_places, true, dist);

            if(status_callback)
                status_callback(std::string("fail_{t")+std::to_string(transitions[0].first)+"}", std::string("nassert"));
            return;
        }

        std::vector<double> probabilities;
        std::vector<std::vector<std::size_t>> all_pre_places;
        std::string transitions_string;

        for(auto [i_transition, probability] : transitions) {
            const auto transition = cn.transitions[i_transition];

            probabilities.push_back(probability);
            all_pre_places.push_back(transition.pre);
            if(status_callback) {
                std::stringstream ss;
                ss << "t" << i_transition << ": " << probability << ", ";
                transitions_string.append(ss.str());
            }
        }

        const auto& pre_places = all_pre_places;
        failp_op(pre_places, probabilities, dist);
        if(status_callback)
            status_callback(std::string("failp_{")+ transitions_string +"}", std::string("failp"));

        return;
    }

    std::vector<double> probabilities;
	std::vector<std::vector<std::size_t>> all_pre_places;
	std::vector<std::vector<std::size_t>> all_post_places;
	std::string transitions_string;

	for(auto [i_transition, probability] : transitions) {
		const auto transition = cn.transitions[i_transition];

		probabilities.push_back(probability);
		all_pre_places.push_back(transition.pre);
		all_post_places.push_back(transition.post);
		if(status_callback) {
			std::stringstream ss;
			ss << "t" << i_transition << ": " << probability << ", ";
			transitions_string.append(ss.str());
		}
	}

	const auto& pre_places = all_pre_places;
	const auto& post_places = all_post_places;
	successp_op(pre_places, post_places, probabilities, dist);
	if(status_callback)
		status_callback(std::string("successp_{")+ transitions_string +"}", std::string("successp"));

	// update marking in CNU
    const auto &transition = cn.transitions[chosen_transition];
    for (const auto &p : transition.pre)
        cn.m.at(p) = 0;
    for (const auto &p : transition.post)
        cn.m.at(p) = 1;
}

void fire_with_probabilityStoch_on_gbn (CN& cn, GBN& gbn, std::vector<std::pair<std::size_t, double>> transitions, std::size_t chosen_transition, std::function<void(std::string,std::string)> status_callback) {
    if(transitions.empty()) {
        std::cout << "stochastic firing: No transitions have been chosen." << std::endl;
        return;
    }

    std::vector<std::pair<std::size_t, double>> valid_transitions;
    for(auto i_transition : transitions) {
        const auto& transition = cn.transitions[i_transition.first];
        if(check_pre_condition(transition, cn.m)) valid_transitions.push_back(i_transition);
    }

    if(valid_transitions.empty()) {

        if(transitions.size() == 1) {
            const auto& transition = cn.transitions[transitions[0].first];
            const auto& pre_places = transition.pre;

            nassert_op(pre_places, true, gbn);
            check_gbn_integrity(gbn);
            if(status_callback)
                status_callback(std::string("fail_{t")+std::to_string(transitions[0].first)+"}", std::string("nassert_pre"));

            return;
        }

        std::vector<double> probabilities;
        std::vector<std::vector<std::size_t>> all_pre_places;
        std::string transitions_string = "";

        for(auto [i_transition, probability] : transitions) {
            const auto transition = cn.transitions[i_transition];

            probabilities.push_back(probability);
            all_pre_places.push_back(transition.pre);
            if(status_callback) {
                std::stringstream ss;
                ss << "t" << i_transition << ": " << probability << ", ";
                transitions_string.append(ss.str());
            }
        }

        auto& pre_places = all_pre_places;
        failp_op(pre_places, probabilities, gbn);
        check_gbn_integrity(gbn);
        if(status_callback) status_callback(std::string("failStoch{")+ transitions_string +"}", std::string("failStoch"));

        return;
    }

    if(transitions.size() == 1) {
        const auto& transition = cn.transitions[transitions[0].first];
        const auto& pre_places = transition.pre;
        const auto& post_places = transition.post;

        assert_op(pre_places,1, gbn);
        check_gbn_integrity(gbn);
        if(status_callback)
            status_callback(std::string("success_{t")+std::to_string(transitions[0].first)+"}", std::string("assert_pre"));

        set_op(pre_places,0,gbn);
        check_gbn_integrity(gbn);
        if(status_callback)
            status_callback(std::string("success_{t")+std::to_string(transitions[0].first)+"}", std::string("set_pre"));

        set_op(post_places,1,gbn);
        check_gbn_integrity(gbn);
        if(status_callback)
            status_callback(std::string("success_{t")+std::to_string(transitions[0].first)+"}", std::string("set_post"));

        if(valid_transitions.empty())
            return;

        // update marking in CNU
        for(const auto& p : transition.pre)
            cn.m.at(p) = 0;
        for(const auto& p : transition.post)
            cn.m.at(p) = 1;
        return;
    }

    std::vector<double> probabilities;
    std::vector<std::vector<std::size_t>> all_pre_places;
    std::vector<std::vector<std::size_t>> all_post_places;
    std::string transitions_string = "";

    for(auto [i_transition, probability] : transitions) {
        const auto transition = cn.transitions[i_transition];

        probabilities.push_back(probability);
        all_pre_places.push_back(transition.pre);
        all_post_places.push_back(transition.post);
        if(status_callback) {
            std::stringstream ss;
            ss << "t" << i_transition << ": " << probability << ", ";
            transitions_string.append(ss.str());
        }
    }

    auto& pre_places = all_pre_places;
    auto& post_places = all_post_places;
    successStoch_op(pre_places, post_places, probabilities, gbn);
    check_gbn_integrity(gbn);
    if(status_callback)
        status_callback(std::string("successStoch_{")+ transitions_string +"}", std::string("successStoch"));

    // update marking in CNU
    if(status_callback)
        status_callback(std::string("chosen transition: t")+ std::to_string(chosen_transition), std::string("successp"));

    const auto &transition = cn.transitions[chosen_transition];
    for (const auto &p : transition.pre)
        cn.m.at(p) = 0;
    for (const auto &p : transition.post)
        cn.m.at(p) = 1;
}

void fire_with_probabilityStoch_on_joint_dist(CN& cn, JointDist& dist, std::vector<std::pair<std::size_t, double>> transitions, std::size_t chosen_transition, std::function<void(std::string,std::string)> status_callback) {
    if(transitions.empty()) {
        std::cout << "stochastic firing: No transitions have been chosen." << std::endl;
        return;
    }

    std::vector<std::pair<std::size_t, double>> valid_transitions;
    for(auto i_transition : transitions) {
        const auto& transition = cn.transitions[i_transition.first];
        if(check_pre_condition(transition, cn.m)) valid_transitions.push_back(i_transition);
    }

    if(valid_transitions.empty()) {

        if(transitions.size() == 1) {
            const auto& transition = cn.transitions[transitions[0].first];
            const auto& pre_places = transition.pre;

            nassert_op(pre_places, true, dist);

            if(status_callback)
                status_callback(std::string("fail_{t")+std::to_string(transitions[0].first)+"}", std::string("nassert"));
            return;
        }

        std::vector<double> probabilities;
        std::vector<std::vector<std::size_t>> all_pre_places;
        std::string transitions_string;

        for(auto [i_transition, probability] : transitions) {
            const auto transition = cn.transitions[i_transition];

            probabilities.push_back(probability);
            all_pre_places.push_back(transition.pre);
            if(status_callback) {
                std::stringstream ss;
                ss << "t" << i_transition << ": " << probability << ", ";
                transitions_string.append(ss.str());
            }
        }

        const auto& pre_places = all_pre_places;
        failp_op(pre_places, probabilities, dist);
        if(status_callback)
            status_callback(std::string("failStoch_{")+ transitions_string +"}", std::string("failStoch"));

        return;
    }

    std::vector<double> probabilities;
    std::vector<std::vector<std::size_t>> all_pre_places;
    std::vector<std::vector<std::size_t>> all_post_places;
    std::string transitions_string;

    for(auto [i_transition, probability] : transitions) {
        const auto transition = cn.transitions[i_transition];

        probabilities.push_back(probability);
        all_pre_places.push_back(transition.pre);
        all_post_places.push_back(transition.post);
        if(status_callback) {
            std::stringstream ss;
            ss << "t" << i_transition << ": " << probability << ", ";
            transitions_string.append(ss.str());
        }
    }

    const auto& pre_places = all_pre_places;
    const auto& post_places = all_post_places;
    successStoch_op(pre_places, post_places, probabilities, dist);
    if(status_callback)
        status_callback(std::string("successStoch_{")+ transitions_string +"}", std::string("successStoch"));

    if(status_callback)
        status_callback(std::string("chosen transition: t")+ std::to_string(chosen_transition), std::string("successp"));

    const auto &transition = cn.transitions[chosen_transition];
    for (const auto &p : transition.pre)
        cn.m.at(p) = 0;
    for (const auto &p : transition.post)
        cn.m.at(p) = 1;
}