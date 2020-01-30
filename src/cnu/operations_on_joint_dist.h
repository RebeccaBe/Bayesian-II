#pragma once

#include "../joint_dist/joint_dist.h"
#include <vector>

void set_op(const std::vector<std::size_t> places, bool b, JointDist& dist);
void setp_op(const std::size_t place, bool b, JointDist& dist, double probability);
void assert_op(const std::vector<std::size_t> places, bool b, JointDist& dist);
void nassert_op(const std::vector<std::size_t> places, bool b, JointDist& dist);
void successp_op(const std::vector<std::vector<std::size_t>> pre_places, const std::vector<std::vector<std::size_t>> post_places,
                const std::vector<double> probabilities, JointDist& dist);
void successStoch_op(const std::vector<std::vector<std::size_t>> pre_places, const std::vector<std::vector<std::size_t>> post_places,
                 const std::vector<double> probabilities, JointDist& dist);

void assert_non_norm_op(const std::vector<std::size_t> places, bool b, JointDist& dist);
void normalize_op(JointDist& dist);
