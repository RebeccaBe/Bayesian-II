#pragma once

#include "cn.h"

#include <random>

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
);
void update_transitions(CN& cn,std::vector<std::size_t>& pre_list);
void update_container(std::unordered_map<std::size_t, std::vector<std::size_t>>& container, std::vector<std::size_t> pool);
void update_post_list(std::vector<size_t> pre_list, std::vector<size_t>& places_for_post);
std::unordered_map<std::size_t, std::vector<std::size_t>> create_container(std::size_t n_places);
void eliminate_loops(CN& cn);