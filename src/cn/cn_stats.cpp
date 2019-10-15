#include "cn_stats.h"
#include "cn_operations.h"

std::size_t get_n_valid_transitions(CN cn) {
    int n_valid_transitions = 0;
    for(auto transition : cn.transitions) {
        if(check_pre_condition(transition,cn.m))
            n_valid_transitions++;
    }
    return n_valid_transitions;
}

int get_n_transitions(CN cn) {
    return cn.transitions.size();
}

int get_n_markings(CN cn) {
    int n_markings = 0;
    for(auto marking : cn.m) {
        if(marking)
            n_markings++;
    }
    return n_markings;
}

std::vector<std::size_t> get_distribution_pre(CN cn){
    std::vector<std::size_t> distribution_pre(cn.n, 0);
    for(auto transition : cn.transitions) {
        for(auto element : transition.pre) {
            distribution_pre[element]++;
        }
    }
    return distribution_pre;
}

std::vector<std::size_t> get_distribution_post(CN cn){
    std::vector<std::size_t> distribution_post(cn.n, 0);
    for(auto transition : cn.transitions) {
        for(auto element : transition.post)
            distribution_post[element]++;
    }
    return distribution_post;
}

double get_average_pre(CN cn) {
    std::size_t n = 0;
    auto distribution = get_distribution_pre(cn);
    for(std::size_t i : distribution)
        n += i;
    double result = (n * 1.0) / (cn.n * 1.0);
    return result;
}

double get_average_post(CN cn) {
    std::size_t n = 0;
    auto distribution = get_distribution_post(cn);
    for(std::size_t i : distribution)
        n += i;
    double result = (n * 1.0) / (cn.n * 1.0);
    return result;
}

std::vector<std::size_t>get_valid_transitions(CN cn) {
    std::vector<std::size_t> valid_transitions;
    for(std::size_t i_transition = 0; i_transition < cn.transitions.size(); i_transition++) {
        if(check_pre_condition(cn.transitions[i_transition],cn.m))
            valid_transitions.push_back(i_transition);
    }
    return valid_transitions;
}