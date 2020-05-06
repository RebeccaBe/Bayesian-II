#include "randomized_generation_fc.h"

#include <algorithm>
#include <iostream>

CN randomize_fc_cn(
	std::size_t n_places,
	std::size_t n_transitions,
	std::size_t n_min_tokens,
	std::size_t n_max_tokens,
	std::size_t n_min_pre_places,
	std::size_t n_max_pre_places,
	std::size_t n_min_post_places,
	std::size_t n_max_post_places,
	std::mt19937& mt
)
{
	if(n_max_tokens > n_places)
		throw std::logic_error("randomize_cn: n_max_tokens > n_places");

	CN cn;
	cn.n = n_places;

	auto container = create_container(n_places);

	std::vector<std::size_t> places(n_places);
	std::iota(places.begin(), places.end(), 0);

    std::vector<std::size_t> places_for_pre(n_places);
    std::vector<std::size_t> places_for_post(n_places);

    std::uniform_int_distribution<std::size_t> rand_tokens(n_min_tokens, n_max_tokens);
    std::uniform_int_distribution<std::size_t> rand_pre_places(n_min_pre_places, n_max_pre_places);
    std::uniform_int_distribution<std::size_t> rand_post_places(n_min_post_places, n_max_post_places);

	for(std::size_t i_transition = 0; i_transition < n_transitions; i_transition++)
	{
		Transition t;

		auto n_pre_places = rand_pre_places(mt);
		auto n_post_places = rand_post_places(mt);

        std::iota(places_for_post.begin(), places_for_post.end(), 0);

		while(n_pre_places > t.pre.size() && n_pre_places>0) {
		    for(auto place : container) {
		        if(place.second.size() <= (n_pre_places - t.pre.size())) {
                    auto res1 = std::find(std::begin(places_for_pre), std::end(places_for_pre), place.first);
                    auto res2 = std::find(std::begin(t.pre), std::end(t.pre), place.first);
                    if(res1 == std::end(places_for_pre) && res2 == std::end(t.pre))
                        places_for_pre.push_back(place.first);
		        }
		    }

		    if(places_for_pre.empty() && t.pre.size() == 0) {
		        n_pre_places++;
		        continue;
		    } else if (places_for_pre.empty() && t.pre.size() > 0) {
		        n_pre_places--;
		        continue;
		    }
            std::shuffle(places_for_pre.begin(), places_for_pre.end(), mt);

		    std::vector<std::size_t> pool = container[places_for_pre.front()]; //pool of specific place
            for(auto p : pool) {
                t.pre.push_back(p);
            }
            update_transitions(cn, t.pre);
            update_container(container, t.pre);

		    places_for_pre.clear();
		}
        update_post_list(t.pre, places_for_post);

        std::shuffle(places_for_post.begin(), places_for_post.end(),mt);
		t.post = std::vector<std::size_t>(places_for_post.begin(), places_for_post.begin()+n_post_places);
        places_for_post.clear();
        places_for_post.resize(n_places);

		cn.transitions.push_back(t);
	}

	std::shuffle(places.begin(), places.end(),mt);
	auto n_tokens = rand_tokens(mt);
	cn.m = Marking(n_places, false);
	for(std::size_t i = 0; i < std::min(n_tokens, n_places); i++)
		cn.m.at(places.at(i)) = true;

	return cn;
}

void update_transitions(CN& cn, std::vector<std::size_t>& pre_list) {
    for (auto& transition : cn.transitions) {
        for (auto list_element : pre_list) {
            auto res = std::find(std::begin(transition.pre), std::end(transition.pre), list_element);
            if(res != std::end(transition.pre))
                transition.pre=pre_list;
        }
    }
}

void update_container(std::unordered_map<std::size_t, std::vector<std::size_t>>& container, std::vector<std::size_t> pool) {
    for(auto place : pool) { //every place out of pool should have the whole pool associated with it in the KV
        container[place] = pool;
    }
}

void update_post_list(std::vector<size_t> pre_list, std::vector<size_t>& places_for_post) {
    for(auto pre_place : pre_list) {
        auto res = std::find(places_for_post.begin(), places_for_post.end(), pre_place);
        if(res != std::end(places_for_post)) {
            places_for_post.erase(res);
        }
    }
}

std::unordered_map<std::size_t, std::vector<std::size_t>> create_container(std::size_t n_places) {
    std::unordered_map<std::size_t, std::vector<std::size_t>> container;
    for(std::size_t i = 0; i < n_places; i++) {
        container[i] = {i};
    }
    return container;
}

void eliminate_loops(CN& cn) {
    for(auto& t : cn.transitions) {
        for (auto el : t.pre) {
            auto res = std::find(t.post.begin(), t.post.end(), el);
            if (res != std::end(t.post))
                t.post.erase(res);
        }
    }
}