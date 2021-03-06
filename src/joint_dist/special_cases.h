#pragma once

#include "joint_dist.h"

JointDist build_uniform_joint_dist(std::size_t n);

JointDist calculate_marginals (std::size_t, JointDist dist);