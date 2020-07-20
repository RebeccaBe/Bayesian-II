#include "special_cases.h"

#include <cmath>

JointDist build_uniform_joint_dist(std::size_t n)
{
	JointDist dist;

	std::vector<bool> marking(n,false);

	double val = 1.0 / std::pow(2,n);

	dist.insert({ marking, val });

	int i = 0;
	while(static_cast<std::size_t>(i) < n)
	{
		if(!marking.at(i))
		{
			marking.at(i) = true;
			while(i > 0)
				marking.at(--i) = false;
			i = 0;

			dist.insert({ marking, val });
		}
		else 
		{
			i++;
		}
	}

	return dist;
}

JointDist calculate_marginals (std::size_t place, JointDist dist) {
    JointDist dist_marg;
    double val_zero = 0;
    double val_one = 0;

    for(auto& [marking,p] : dist) {
        if(marking[place] == 0)
            val_zero += p;
        else
            val_one += p;
    }

    dist_marg.insert({std::vector<bool> (1,false), val_zero});
    dist_marg.insert({std::vector<bool> (1,true), val_one});

    return dist_marg;
}
