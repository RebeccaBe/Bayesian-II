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

    /*auto i_valid_transitions = get_valid_transitions(cn);
    std::shuffle(i_valid_transitions.begin(), i_valid_transitions.end(),mt);
    std::size_t valid_transition = i_valid_transitions.at(0); //guarantee one valid transition

    auto rest_transitions = i_transitions;
    rest_transitions.erase(rest_transitions.begin()+valid_transition);
    std::shuffle(rest_transitions.begin(), rest_transitions.end(),mt);
    auto chosen_transitions = std::vector<std::size_t>(rest_transitions.begin(), rest_transitions.begin()+n_transitions-1);
    chosen_transitions.push_back(valid_transition);*/
    std::shuffle(i_transitions.begin(), i_transitions.end(),mt);
    auto chosen_transitions = std::vector<std::size_t>(i_transitions.begin(), i_transitions.begin() + n_transitions);

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

std::vector<std::pair<std::size_t, double>> RandomTransitionHelper::next_from_bubbles(std::mt19937& mt) {
    if(type != PROBABILITY)
        throw std::logic_error("Wrong transition type, cannot perform probability operation.");

    /*std::vector<std::vector<std::pair<std::size_t, double>>> bubbles_w_valid_transitions;
    for(auto i_bubble : transition_bubbles) {
        for(auto [i_transition, p] : i_bubble) {
            if(check_pre_condition(cn.transitions[i_transition], cn.m)) {
                bubbles_w_valid_transitions.push_back(i_bubble);
                break; //for efficiency
            }
        }
    }*/

    std::shuffle(transition_bubbles.begin(), transition_bubbles.end(), mt);
    return  transition_bubbles[0];
}

// This is only needed for correct testing. If we would not choose the transition beforehand, maybe different ones would be chosen for
// the same test and therefore, the resulting joint distributions could not be compared to each other.
std::size_t RandomTransitionHelper::choose_transition (CN& cn, std::vector<std::pair<std::size_t, double>> transitions_w_probabilities) {

    std::vector<std::pair<std::size_t, double>> valid_transitions;
    for(auto i_transition : transitions_w_probabilities) {
        const auto& transition = cn.transitions[i_transition.first];
        if(check_pre_condition(transition, cn.m)) valid_transitions.push_back(i_transition);
    }

    if(valid_transitions.empty())
        return cn.transitions.size(); //No valid transition

    std::random_device rd;
    std::mt19937 mt(rd());

    double normalization = 0;
    for (auto t : valid_transitions) normalization += t.second;

    std::uniform_real_distribution<double> rand_transition_0_1(0, 1);
    auto chosen_transition_p = rand_transition_0_1(mt);

    double sum = 0;
    std::size_t index = 0;
    while (sum <= chosen_transition_p) {
        sum += (valid_transitions[index].second / normalization);
        index++;
    }

    return valid_transitions[index-1].first;
}

std::vector<std::vector<std::pair<std::size_t, double>>> RandomTransitionHelper::make_transitions_w_probabilities(std::mt19937& mt) { //new class for this?
    std::vector<std::vector<std::pair<std::size_t, double>>> transitions_w_probabilities_bubbles;

    if(type != PROBABILITY) return transitions_w_probabilities_bubbles;

    for(std::size_t i_bubble = 0; i_bubble < cn.transitions.size()/3; i_bubble++) {
        std::vector<std::pair<std::size_t, double>> i_bubble_vec;

        std::vector<std::size_t> chosen_transitions;
        std::vector<double> chosen_probabilities;

        std::uniform_int_distribution<std::size_t> rand_n_transitions(1, max_trans);
        auto n_transitions = rand_n_transitions(mt);
        if(n_transitions > cn.transitions.size()) n_transitions = cn.transitions.size();

        std::shuffle(i_transitions.begin(), i_transitions.end(),mt);
        chosen_transitions = std::vector<std::size_t>(i_transitions.begin(), i_transitions.begin() + n_transitions);

        for(std::size_t i = 0; i < n_transitions-1; i++)
            chosen_probabilities.push_back(uniform_0_1(mt));
        std::sort(chosen_probabilities.begin(), chosen_probabilities.end());

        double probability_sum = 0;
        for(std::size_t i = 0; i < n_transitions-1; i++) {
            i_bubble_vec.emplace_back(chosen_transitions[i],(chosen_probabilities[i]-probability_sum));
            probability_sum = chosen_probabilities[i];
        }
        i_bubble_vec.emplace_back(chosen_transitions.back(),1-probability_sum);

        transitions_w_probabilities_bubbles.push_back(i_bubble_vec);
    }

    return transitions_w_probabilities_bubbles;
}