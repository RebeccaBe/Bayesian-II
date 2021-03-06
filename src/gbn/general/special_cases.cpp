#include "special_cases.h"

GBN build_uniform_independent_obn(std::size_t n_vertices)
{
	GBN gbn(0,n_vertices,n_vertices);
	auto& g = gbn.graph;

	for(auto v : ::inside_vertices(gbn))
	{
		auto p_m = std::make_shared<DynamicMatrix>(0,1);
		p_m->set(0,0,0.5);
		p_m->set(1,0,0.5);
		put(vertex_matrix, g, v, p_m);

		auto e_new = boost::add_edge(v,v+n_vertices,g).first;
		put(edge_position, g, e_new, std::make_pair(0,0));
        put(edge_equivalence_class, g, e_new, max_equivalence_counter_and_increase(gbn));
	}

	return gbn;
}

GBN build_independent_obn(std::vector<std::pair<double,double>> rv_dists)
{
	GBN gbn(0,rv_dists.size(),rv_dists.size());
	auto& g = gbn.graph;

	for(auto v : ::inside_vertices(gbn))
	{
		auto p_m = std::make_shared<DynamicMatrix>(0,1);
		p_m->set(0,0,rv_dists[v].first);
		p_m->set(1,0,rv_dists[v].second);
		put(vertex_matrix, g, v, p_m);

		auto e_new = boost::add_edge(v,v+rv_dists.size(),g).first;
		put(edge_position, g, e_new, std::make_pair(0,0));
        put(edge_equivalence_class, g, e_new, max_equivalence_counter_and_increase(gbn));
	}

	return gbn;
}

GBN build_obn_without_knowledge(std::size_t n_vertices)
{
	GBN gbn(n_vertices, n_vertices, 0);
	auto& g = gbn.graph;

	for(auto v : ::input_vertices(gbn))
	{
		auto e_new = boost::add_edge(v,v+n_vertices,g).first;
		put(edge_position, g, e_new, std::make_pair(0,0));
        put(edge_equivalence_class, g, e_new, max_equivalence_counter_and_increase(gbn));
	}

	return gbn;
}
