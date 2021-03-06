#include "cn_parameters.h"

CN_Parameters read_cn_params(const cxxopts::ParseResult& params)
{
	CN_Parameters rtn;

	rtn.N_MIN_PLACES = params["n-min-places"].as<std::size_t>();
	rtn.N_MAX_PLACES = params["n-max-places"].as<std::size_t>();
	rtn.PERCENT_TRANSITIONS = params["percent-transitions"].as<std::size_t>();
	rtn.N_RUNS = params["n-runs"].as<std::size_t>();
	rtn.N_TRANSITIONS_PER_RUN = params["n-transitions-per-run"].as<std::size_t>();
	rtn.PERCENT_MIN_TOKENS = params["percent-min-tokens"].as<double>();
	rtn.PERCENT_MAX_TOKENS = params["percent-max-tokens"].as<double>();

	rtn.N_MIN_PRE_PLACES = params["n-min-pre-places"].as<std::size_t>();
	rtn.N_MAX_PRE_PLACES = params["n-max-pre-places"].as<std::size_t>();
	rtn.N_MIN_POST_PLACES = params["n-min-post-places"].as<std::size_t>();
	rtn.N_MAX_POST_PLACES = params["n-max-post-places"].as<std::size_t>();

	rtn.RAND_TRANSITION_TYPE = (params["rand-transition-type"].as<std::string>() == "forced") ?
	        RandomTransitionHelper::FORCED : (params["rand-transition-type"].as<std::string>() == "probability") ?
	        RandomTransitionHelper::PROBABILITY : RandomTransitionHelper::UNIFORM;
	rtn.P_SUCCESS = params["p-success"].as<double>();

	rtn.FREE_CHOICE = params["free-choice"].as<bool>();
    rtn.EVALUATION_TYPE = (params["evaluation_type"].as<std::string>() == "degree") ?
            DEGREE : (params["evaluation_type"].as<std::string>() == "fillIn") ?
            FILLIN : DEFAULT;
	//rtn.COMPARISON_VALUE = params["comparison-value"].as<std::size_t>();

	rtn.N_MAX_TRANSITIONS_PER_OP = params["n-max-transitions-per-op"].as<std::size_t>();

	return rtn;
}

void write_cn_params(std::ostream& ostr, const CN_Parameters& cn_params)
{
	ostr << "N_MIN_PLACES: " << cn_params.N_MIN_PLACES << std::endl;
	ostr << "N_MAX_PLACES: " << cn_params.N_MAX_PLACES << std::endl;
	ostr << "N_PERCENT_TRANSITIONS: " << cn_params.PERCENT_TRANSITIONS << std::endl;
	ostr << "N_RUNS: " << cn_params.N_RUNS << std::endl;
	ostr << "N_TRANSITIONS_PER_RUN: " << cn_params.N_TRANSITIONS_PER_RUN << std::endl;
	ostr << "PERCENT_MIN_TOKENS: " << cn_params.PERCENT_MIN_TOKENS << std::endl;
	ostr << "PERCENT_MAX_TOKENS: " << cn_params.PERCENT_MAX_TOKENS << std::endl;

	ostr << "N_MIN_PRE_PLACES: " << cn_params.N_MIN_PRE_PLACES << std::endl;
	ostr << "N_MAX_PRE_PLACES: " << cn_params.N_MAX_PRE_PLACES << std::endl;
	ostr << "N_MIN_POST_PLACES: " << cn_params.N_MIN_POST_PLACES << std::endl;
	ostr << "N_MAX_POST_PLACES: " << cn_params.N_MAX_POST_PLACES << std::endl;

	ostr << "RAND_TRANSITION_TYPE: " << (cn_params.RAND_TRANSITION_TYPE == RandomTransitionHelper::FORCED ?
	    "FORCED" : (cn_params.RAND_TRANSITION_TYPE == RandomTransitionHelper::PROBABILITY) ?
	    "PROBABILITY" : "UNIFORM") << std::endl;
	ostr << "P_SUCCESS: " << cn_params.P_SUCCESS << std::endl;

	ostr << "FREE_CHOICE: " << cn_params.FREE_CHOICE << std::endl;
    ostr << "EVALUATION_TYPE: " << (cn_params.EVALUATION_TYPE == DEGREE ?
        "DEGREE" : (cn_params.EVALUATION_TYPE == FILLIN) ?
        "FILLIN" : "DEFAULT") << std::endl;
	//ostr << "COMPARISON_VALUE: " << cn_params.COMPARISON_VALUE << std::endl;

	ostr << "N_MAX_TRANSITIONS_PER_OP: " << cn_params.N_MAX_TRANSITIONS_PER_OP << std::endl;
}

