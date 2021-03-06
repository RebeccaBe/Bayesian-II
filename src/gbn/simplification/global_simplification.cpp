#include "simplification.h"

#include <iostream>
#include "../matrix/matrix_io.h"
#include "../modification/vertex_add_remove.h"
#include "../modification/merging.h"
#include "../modification/splitting.h"
#include "../general/path_closing.h"
#include "../general/gbn_io.h"
#include "../general/check.h"
#include "../../helpers.hpp"

bool eliminate_stochastic_vertex_without_outputs(GBN& gbn, Vertex v, std::string& op)
{
	op = "eliminate_stochastic_vertex_without_outputs";
	auto& g = gbn.graph;

	if(matrix(v,g)->type == TERMINATOR || !matrix(v,g)->is_stochastic)
		return false;

	bool connected_to_output = false;
	for(auto e : boost::make_iterator_range(boost::out_edges(v,g)))
	{
		if(type(boost::target(e,g),g) == OUTPUT)
		{
			connected_to_output = true;
			break;
		}
	}

	bool has_sucessors = false;
	if(!connected_to_output) {
		std::vector<Vertex> successors;
		successors.push_back(v);
		for(auto e : boost::make_iterator_range(boost::out_edges(v,g)))
			successors.push_back(boost::target(e,g));

		has_sucessors = successors.size() > 1;
		if(has_sucessors) {
            //TODO: check predecessors of successors for efficiency purposes


            successors = path_closing(gbn, successors);

			auto v_new = merge_vertices(gbn, successors);

			//TODO?: think about a kind of case distinction, when the vertex should be split and when not...?
			//recursively_split_vertex(gbn, v_new);
		}
	}

	return !connected_to_output && has_sucessors;
}

bool switch_substoch_to_front(GBN& gbn, Vertex v, std::string& op)
{
	op = "switch_substoch_to_front";
	auto& g = gbn.graph;

	if(matrix(v,g)->is_stochastic)
		return false;

	// check if is already in front
	bool only_inputs_as_predecessors = true;
	for(auto e : boost::make_iterator_range(boost::in_edges(v,g)))
		if(type(boost::source(e,g), g) != INPUT)
		{
			only_inputs_as_predecessors = false;
			break;
		}

	if(only_inputs_as_predecessors)
		return false;

	std::vector<Vertex> precessors;
	precessors.push_back(v);
	for(auto e : boost::make_iterator_range(boost::in_edges(v,g)))
	{
		auto u = boost::source(e,g);
		if(type(u,g) == NODE)
			precessors.push_back(u);
	}

	precessors = path_closing(gbn, precessors);
	auto v_new = merge_vertices(gbn, precessors);

	recursively_split_vertex(gbn, v_new);

	return true;
}

bool normalize_substoch_front_vertices_without_inputs(GBN& gbn, Vertex v, std::string& op)
{
	op = "normalize_substoch_front_vertices_without_inputs";
	auto& g = gbn.graph;

	if(matrix(v,g)->is_stochastic || boost::in_degree(v,g) > 0 || matrix(v,g)->type == ZERO)
		return false;

	auto m = matrix(v,g);

	// special cases: if only one entry not equal zero, renormalize to 1_0 resp. 1_1
	if(m->m <= 1 && prob_is_zero(m->get(BitVec(0),BitVec(0))) && !prob_is_zero(m->get(BitVec(1),BitVec(0))))
	{
		auto p_m_new = std::make_shared<OneBMatrix>(true);
		put(vertex_matrix, g, v, p_m_new);
		return true;
	}

	if(m->m <= 1 && !prob_is_zero(m->get(BitVec(0),BitVec(0))) && prob_is_zero(m->get(BitVec(1),BitVec(0))))
	{
		auto p_m_new = std::make_shared<OneBMatrix>(false);
		put(vertex_matrix, g, v, p_m_new);
		return true;
	}

	/*double sum = m->get(BitVec(0),BitVec(0)) + m->get(BitVec(1),BitVec(0));

	m->set(BitVec(0),BitVec(0), m->get(BitVec(0),BitVec(0))/sum);
	m->set(BitVec(1),BitVec(0), m->get(BitVec(1),BitVec(0))/sum);*/

	double sum = 0;
	for(std::size_t assignmentX = 0; assignmentX <= pow(2, (m->m)); assignmentX++)
		sum += m->get(assignmentX, BitVec(0));

	for(std::size_t assignmentX = 0; assignmentX <= pow(2, (m->m)); assignmentX++)
        m->set(assignmentX, BitVec(0), m->get(assignmentX, BitVec(0)) / sum);

	m->is_stochastic = true;

	return true;
}
