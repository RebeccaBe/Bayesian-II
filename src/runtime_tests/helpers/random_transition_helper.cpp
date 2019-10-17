#include "random_transition_helper.h"

RandomTransitionHelper::RandomTransitionHelper(const CN& cn, Type type, double p_success, std::size_t max_trans)
	:cn(cn), 
	i_transitions(cn.transitions.size(), 0),
	type(type), 
	p_success(p_success),
	max_trans(max_trans),
	rand_transition(0,cn.transitions.size()-1),
	uniform_0_1(0,1),
     n_valid_transitions(get_n_valid_transitions(cn))
{
	std::iota(i_transitions.begin(), i_transitions.end(), 0);
}

std::size_t RandomTransitionHelper::next(std::mt19937& mt) 
{
	switch(type) {
		case UNIFORM:
			return rand_transition(mt);

		case FORCED: 
			{
				bool force_success = p_success > uniform_0_1(mt);

				if(force_success) {
					std::shuffle(i_transitions.begin(), i_transitions.end(), mt);
					for(auto i_transition : i_transitions)
					{
						if(check_pre_condition(cn.transitions[i_transition],cn.m))
							return i_transition;
					}
					return i_transitions.front();
				}
				else
					return rand_transition(mt);
			}
		default:
			throw std::logic_error("Undefined case.");
	}
}

std::vector<std::pair<std::size_t, double>> RandomTransitionHelper::next_p(std::mt19937& mt) {
	if(type != PROBABILITY)
		throw std::logic_error("Wrong transition type, cannot perform probability operation.");

    std::vector<std::pair<std::size_t, double>> transitions_w_probabilities;

    n_valid_transitions=get_n_valid_transitions(cn);

    if(n_valid_transitions < 1)
        return transitions_w_probabilities;

	std::uniform_int_distribution<std::size_t> rand_n_transitions(1, max_trans);
	auto n_transitions = rand_n_transitions(mt);
    if(n_transitions > cn.transitions.size()) n_transitions=cn.transitions.size();

    auto i_valid_transitions = get_valid_transitions(cn);
    std::shuffle(i_valid_transitions.begin(), i_valid_transitions.end(),mt);
    std::size_t valid_transition = i_valid_transitions.at(0); //guarantee one valid transition

    auto rest_transitions = i_transitions;
    rest_transitions.erase(rest_transitions.begin()+valid_transition);
    std::shuffle(rest_transitions.begin(), rest_transitions.end(),mt);
    auto chosen_transitions = std::vector<std::size_t>(rest_transitions.begin(), rest_transitions.begin()+n_transitions-1);
    chosen_transitions.push_back(valid_transition);

	std::vector<double> probabilities;
    for(std::size_t i = 0; i < n_transitions-1; i++)
        probabilities.push_back(uniform_0_1(mt));
    std::sort(probabilities.begin(), probabilities.end());

	double probability = 0;
	for(std::size_t i = 0; i < n_transitions-1; i++) {
        transitions_w_probabilities.emplace_back(chosen_transitions[i],(probabilities[i]-probability));
        probability = probabilities[i];
	}
    transitions_w_probabilities.emplace_back(chosen_transitions.back(),1-probability);

	return transitions_w_probabilities;
}