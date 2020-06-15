#include "simplification.h"
#include "local_simplification.h"

#include <iostream>
#include <map>
#include "../matrix/matrix_io.h"
#include "../modification/vertex_add_remove.h"
#include "../modification/merging.h"
#include "../modification/splitting.h"
#include "../general/path_closing.h"
#include "../general/gbn_io.h"
#include "../../helpers.hpp"

#include <fstream>
#include <unistd.h>

bool check_and_apply_F1(GBN& gbn, Vertex v, std::string& op)
{
	op = "check_and_apply_F1";
	auto& g = gbn.graph;

	if(type(v,g) != NODE)
		return false;

	if(matrix(v,g)->type == ONE_B)
	{
		if(boost::out_degree(v,g) > 1)
		{
			auto& m = dynamic_cast<OneBMatrix&>(*matrix(v,g));

			std::vector<Port> sucessor_ports;	
			for(auto e : boost::make_iterator_range(boost::out_edges(v,g)))
				sucessor_ports.push_back({ boost::target(e,g), port_to(e,g) });

			for(auto p : sucessor_ports)
			{
				auto v_new = add_vertex(gbn, std::make_shared<OneBMatrix>(m.b),std::string("1_") + std::to_string(m.b));
				auto e = boost::add_edge(v_new, p.first, g).first;
				auto eq_class = equivalence_class(get_edge(std::pair<Vertex, std::size_t>{v, 0}, std::pair<Vertex, std::size_t>{p.first, p.second}, g), g);
				put(edge_position, g, e, std::pair<std::size_t, std::size_t>{ 0, p.second });
                put(edge_equivalence_class, g, e, eq_class);
			}

			remove_vertex(v,gbn);

			return true;
		}
	}

	return false;
}

bool check_and_apply_F2(GBN& gbn, Vertex v, std::string& op)
{
	op = "check_and_apply_F2";
	auto& g = gbn.graph;

	if(type(v,g) != NODE)
		return false;

	if(!matrix(v,g)->is_stochastic)
		return false;

	if(boost::out_degree(v,g) == 0)
		return false;

	for(auto e : boost::make_iterator_range(boost::out_edges(v,g)))
	{
		auto v_suc = boost::target(e,g);
		if(type(v_suc,g) != NODE)
			return false;
		if(matrix(v_suc,g)->type != TERMINATOR)
			return false;
	}

	// if we reached here, then v_pre is a stochastic vertex whose only successor is a terminator

	std::vector<Port> precessor_ports;
	std::map<Port, size_t> eq_classes;
	for(auto e : boost::make_iterator_range(boost::in_edges(v,g))) {
        precessor_ports.push_back({boost::source(e, g), port_from(e, g)});
        eq_classes[{boost::source(e, g), port_from(e, g)}] = equivalence_class(e, g);
    }

	for(auto e : boost::make_iterator_range(boost::out_edges(v,g)))
	{
		auto v_suc = boost::target(e,g);
		remove_vertex(v_suc,gbn);
	}

	for(auto& p : precessor_ports)
	{
		auto v_term = add_vertex(gbn, std::make_shared<TerminatorMatrix>(), "T");
		auto e = boost::add_edge(p.first, v_term, g).first;
		put(edge_position, g, e, std::pair<std::size_t, std::size_t>{ p.second, 0 });
        put(edge_equivalence_class, g, e, eq_classes[p]);
	}

    remove_vertex(v,gbn);

	return true;
}

// TODO: implement the more general version for n -> m
bool check_and_apply_CoUnit(GBN& gbn, Vertex v, std::string& op)
{
	op = "check_and_apply_CoUnit";
	auto& g = gbn.graph;

	if(type(v,g) != NODE || matrix(v,g)->type != TERMINATOR)
		return false;

	auto e_pre = *(boost::in_edges(v,g).first);
	auto v_pre = boost::source(e_pre, g);

	auto term_pre_port = port_from(e_pre, g);
	std::size_t n_same_ports = 0;
	for(auto e : boost::make_iterator_range(boost::out_edges(v_pre,g)))
	{
		if(port_from(e,g) == term_pre_port)
			n_same_ports++;
	}

	if(n_same_ports <= 1)
		return false;
	else
		remove_vertex(v,gbn);

	return true;
}

bool check_and_apply_F3(GBN& gbn, Vertex v, std::string& op)
{
	op = "check_and_apply_F3";
	auto& g = gbn.graph;

	if(type(v,g) != NODE || matrix(v,g)->type != F)
		return false;

	// list of input ports of FMatrix node
	std::map<Port, std::size_t> input_ports_of_F;
	for(auto e : boost::make_iterator_range(boost::in_edges(v,g)))
		input_ports_of_F.insert({{ boost::source(e,g), port_from(e,g) }, port_to(e,g) });

	// list successors of Fmatrix
	std::set<Vertex> successors;
	for(auto e : boost::make_iterator_range(boost::out_edges(v,g)))
		successors.insert(boost::target(e,g));

	// run through all successors 
	bool found_once = false;
	bool found;
	do {
		found = false;
		Port pre_port, post_port;
		Index F_idx;

		for(auto v_suc : successors)
		{
			for(auto [it, end_it] = boost::in_edges(v_suc,g); it != end_it && !found; it++)
			{
				auto e = *it;
				pre_port = { boost::source(e,g), port_from(e,g) };
				post_port = { boost::target(e,g), port_to(e,g) };

				auto map_it = input_ports_of_F.find(pre_port);
				if(map_it != input_ports_of_F.end())
				{
					F_idx = map_it->second;
					found = true;
				}
			}
			if(found)
				break;
		}

		if(found) {
			found_once = true;
			auto eq_class = equivalence_class(get_edge(pre_port, post_port, g), g); // TODO: correct??
			boost::remove_edge_if([&](const Edge& e) {
				return boost::source(e,g) == pre_port.first && boost::target(e,g) == post_port.first && port_from(e,g) == pre_port.second && port_to(e,g) == post_port.second;	
			}, g);
			auto e_new = boost::add_edge(v,post_port.first,g).first;
			put(edge_position, g, e_new, std::pair<std::size_t, std::size_t>{ F_idx, post_port.second });
            put(edge_equivalence_class, g, e_new, eq_class);
		}
	}
	while(found);

	return found_once;
}

bool check_and_apply_F4(GBN& gbn, Vertex v_oneb, std::string& op)
{
	op = "check_and_apply_F4";
	auto& g = gbn.graph;

	if(type(v_oneb,g) != NODE || matrix(v_oneb,g)->type != ONE_B || boost::out_degree(v_oneb, g) != 1)
		return false;

	auto& m = dynamic_cast<OneBMatrix&>(*matrix(v_oneb,g));
	auto b = m.b;

	bool found = false;
	Vertex v_F;
	Index F_port;
	for(auto e : boost::make_iterator_range(boost::out_edges(v_oneb,g)))
	{
		auto v_to = boost::target(e,g);
		if(type(v_to,g) == NODE && matrix(v_to, g)->type == F)
		{
			auto& m_F = dynamic_cast<FMatrix&>(*matrix(v_to, g));
			if(m_F.b == b)
			{
				v_F = v_to;
				F_port = port_to(e,g);
				found = true;
				break;
			}
		} /*else if (type(v_to,g) == NODE && matrix(v_to, g)->type == DIAGONAL) {
			auto& m_DIAG = dynamic_cast<DiagonalMatrix&>(*matrix(v_to, g));
			if((m_DIAG.ones.size() != (pow(2, m_DIAG.k)-1)) ||
				(m_DIAG.get(0,0) == 1 && m_DIAG.get((pow(2, m_DIAG.k)-1),(pow(2, m_DIAG.k)-1)) == 1)) {
				continue;
			}

			bool zero_entry = (m_DIAG.get(0,0) == 0) ? 0 : 1;
			if(zero_entry == b) {
				auto m_F_new = std::make_shared<FMatrix>(m_DIAG.k, zero_entry);
				put(vertex_matrix, gbn.graph, v_to, m_F_new);

				v_F = v_to;
				F_port = port_to(e,g);
				found = true;
				break;
			}
		}*/
	}

	if(!found) 
		return false;

	auto p_m_F = matrix(v_F, g);
	auto& m_F = dynamic_cast<FMatrix&>(*p_m_F);

	std::vector<Port> output_ports;
	std::map<Port, std::size_t> eq_classes;
	for(auto e : boost::make_iterator_range(boost::out_edges(v_F,g)))
		if(port_from(e,g) == F_port) {
            output_ports.push_back({boost::target(e, g), port_to(e, g)});
            eq_classes[{boost::target(e, g), port_to(e, g)}] = equivalence_class(e, g);
        }

	if(m_F.k == 1) 
	{
		remove_vertex(v_oneb,gbn);
		remove_vertex(v_F,gbn);
		auto p_m_new = std::make_shared<ZeroMatrix>(0,1);
		auto v = add_vertex(gbn, p_m_new, "0");
		for(auto p : output_ports)
		{
			auto e = boost::add_edge(v,p.first,g).first;
			put(edge_position,g,e,std::pair<std::size_t,std::size_t>{ 0, p.second });
            put(edge_equivalence_class,g,e,eq_classes[p]);
		}
	}
	else
	{
		auto p_m_F_new = std::make_shared<FMatrix>(m_F.k-1, m_F.b);

		boost::remove_edge_if([&](const Edge& e) {
				if(boost::target(e,g) == v_F && port_to(e,g) == F_port)
				return true;
				if(boost::source(e,g) == v_F && port_from(e,g) == F_port)
				return true;

				return false;
				},g);

		for(auto e : boost::make_iterator_range(boost::in_edges(v_F,g)))
			if(port_to(e,g) >= F_port)
				put(edge_position, g, e, std::pair<std::size_t, std::size_t>{ port_from(e,g), port_to(e,g)-1 });
		for(auto e : boost::make_iterator_range(boost::out_edges(v_F,g)))
			if(port_from(e,g) >= F_port)
				put(edge_position, g, e, std::pair<std::size_t, std::size_t>{ port_from(e,g)-1, port_to(e,g) });

		put(vertex_matrix, g, v_F, p_m_F_new);

		for(auto p : output_ports)
		{
			auto e = boost::add_edge(v_oneb,p.first,g).first;
			put(edge_position,g,e,std::pair<std::size_t,std::size_t>{ 0, p.second });
            put(edge_equivalence_class,g,e,eq_classes[p]);
		}
	}

	return true;
}

bool check_and_apply_F5(GBN& gbn, Vertex v_oneb, std::string& op)
{
	op = "check_and_apply_F5";
	auto& g = gbn.graph;

	if(type(v_oneb,g) != NODE || matrix(v_oneb,g)->type != ONE_B || boost::out_degree(v_oneb, g) != 1)
		return false;

	auto& m_oneb = dynamic_cast<OneBMatrix&>(*matrix(v_oneb,g));

	bool found = false;
	Vertex v_F;
	Index F_port;
	for(auto e : boost::make_iterator_range(boost::out_edges(v_oneb,g)))
	{
		auto v_to = boost::target(e,g);
		if(type(v_to,g) == NODE && matrix(v_to, g)->type == F)
		{
			auto& m_F = dynamic_cast<FMatrix&>(*matrix(v_to, g));
			if(m_F.k < 2)
				continue;
			if(m_F.b == !m_oneb.b)
			{
				v_F = v_to;
				F_port = port_to(e,g);
				found = true;
				break;
			}
		} /*else if (type(v_to,g) == NODE && matrix(v_to, g)->type == DIAGONAL) {
			auto& m_DIAG = dynamic_cast<DiagonalMatrix&>(*matrix(v_to, g));
			if((m_DIAG.ones.size() != (pow(2, m_DIAG.k)-1)) ||
			   (m_DIAG.get(0,0) == 1 && m_DIAG.get((pow(2, m_DIAG.k)-1),(pow(2, m_DIAG.k)-1)) == 1)) {
				continue;
			}

			bool zero_entry = (m_DIAG.get(0,0) == 0) ? 0 : 1;
			if(zero_entry == !m_oneb.b) {
                auto m_F_new = std::make_shared<FMatrix>(m_DIAG.k, zero_entry);
                put(vertex_matrix, gbn.graph, v_to, m_F_new);

                v_F = v_to;
                F_port = port_to(e,g);
                found = true;
                break;
			}
		}*/
	}

	if(!found) 
		return false;

	auto& m_F = dynamic_cast<FMatrix&>(*matrix(v_F, g));

	for(std::size_t i_port = 0; i_port < m_F.k; i_port++)
	{
		if(i_port == F_port)
			continue;

		Vertex v_pre;
		Index port_pre;
		for(auto e : boost::make_iterator_range(boost::in_edges(v_F,g)))
		{
			if(port_to(e,g) == i_port)
			{
				v_pre = boost::source(e,g);
				port_pre = port_from(e,g);
				break;
			}
		}

		for(auto e : boost::make_iterator_range(boost::out_edges(v_F,g)))
		{
			if(port_from(e,g) == i_port)
			{
				auto e_new = boost::add_edge(v_pre, boost::target(e,g), g).first;
				put(edge_position, g, e_new, std::pair<std::size_t, std::size_t>{ port_pre, port_to(e,g) });
                put(edge_equivalence_class, g, e_new, equivalence_class(e, g));
			}
		}
	}
	for(auto e : boost::make_iterator_range(boost::out_edges(v_F,g)))
	{
		if(port_from(e,g) == F_port)
		{
			auto e_new = boost::add_edge(v_oneb, boost::target(e,g), g).first;	
			put(edge_position, g, e_new, std::pair<std::size_t, std::size_t>{ 0, port_to(e,g) });
            put(edge_equivalence_class, g, e_new, equivalence_class(e, g));
		}
	}

	remove_vertex(v_F,gbn);

	return true;
}

bool split_vertex_if_multiple_outputs(GBN& gbn, Vertex v, std::string& op)
{
	op = "split_vertex_if_multiple_outputs";
	auto& g = gbn.graph;

	auto& m = *matrix(v,g);

	if(m.m <= 1 || !m.is_stochastic)
		return false;

	recursively_split_vertex(gbn, v);

	return true;
}

bool simplify_matrix_for_duplicate_inputs(GBN& gbn, Vertex v, std::string& op)
{
	op = "simplify_matrix_for_duplicate_inputs";
	auto& g = gbn.graph;
	auto& m = *matrix(v,g);

	std::map<Port,std::vector<std::size_t>> external_input_to_input_port;
    std::map<Port,std::size_t> eq_classes;
	for(auto e : boost::make_iterator_range(boost::in_edges(v,g)))
	{
		auto ext_input = std::make_pair(boost::source(e,g), port_from(e,g));
		external_input_to_input_port[ext_input].push_back(port_to(e,g));
		eq_classes[ext_input] = equivalence_class(e, g);
	}

	bool found_duplicate = false;
	std::vector<std::vector<std::size_t>> new_to_old_map;
	for(auto t : external_input_to_input_port)
	{
		if(t.second.size() > 1)
			found_duplicate = true;

		new_to_old_map.push_back(t.second);	
	}

	if(!found_duplicate)
		return false;

	new_to_old_map.erase(std::remove_if(new_to_old_map.begin(), new_to_old_map.end(), [](const auto& vec) { return vec.empty(); }), new_to_old_map.end());

	// build new matrix TODO: case of F matrix
	MatrixPtr p_m_new;
	if(m.type == F) {
		auto& m2 = dynamic_cast<FMatrix&>(m);
		p_m_new = std::make_shared<FMatrix>(new_to_old_map.size(), m2.b);
	} else {
		p_m_new = std::make_shared<DynamicMatrix>(new_to_old_map.size(), m.m);
		auto& m_new = *p_m_new;

		unsigned long long i_max_row = 1;
		unsigned long long i_max_col = 1;
		i_max_col = i_max_col << m_new.n;
		i_max_row = i_max_row << m_new.m;
		for(Index i_row = 0; i_row < i_max_row; i_row++)
			for(Index i_col = 0; i_col < i_max_col; i_col++)
			{
				BitVec to(i_row);
				BitVec from_new(i_col);
				BitVec from_old;
				from_old.reset();

				for(std::size_t i = 0; i < new_to_old_map.size(); i++)
					for(auto i_port : new_to_old_map[i])
						if(from_new[i])
							from_old.set(i_port);

				m_new.set(to, from_new, m.get(to,from_old));
			}

		put(vertex_matrix, g, v, p_m_new);
	}

	// rewire vertex
	boost::clear_in_edges(v,g);

	std::size_t i_counter = 0;
	for(const auto [port, i_ports] : external_input_to_input_port)
	{
		auto e = boost::add_edge(port.first, v, g).first;
		put(edge_position, g, e, std::make_pair(port.second,i_counter));
        put(edge_equivalence_class, g, e, eq_classes[port]);
		i_counter++;
	}

	return true;
}

bool merge_F_matrices_to_diagonal_matrix(GBN& gbn, Vertex v, std::string& op) {

	op = "merge_F_matrices_to_diagonal_matrix";

	auto& g = gbn.graph;
	auto& m = *matrix(v,g);

	if(!(m.type == F || m.type == DIAGONAL)) return false;
	else {

	    auto vecs_out = find_recursively_out(gbn, v);
	    auto vecs_in = find_recursively_in(gbn, v);

	    std::set<Vertex> mergeable_vertices = vecs_out;
	    mergeable_vertices.insert(vecs_in.begin(), vecs_in.end());

        bool found_mergeable_vertices = (mergeable_vertices.size() > 1) ? true : false;

		if(found_mergeable_vertices) {
			merge_vertices(gbn, std::vector<Vertex>(mergeable_vertices.begin(), mergeable_vertices.end()), "diagonal");
			return true;
		} else
			return false;
	}
}

std::set<Vertex> find_recursively_out(GBN gbn, Vertex v) { //can there be loops in a gbn?
	auto g = gbn.graph;
	bool path_closed = true;

	std::set<Vertex> vertices;
	vertices.insert(v);

	for(auto e : boost::make_iterator_range(boost::out_edges(v, g))) {

		Vertex v_to = boost::target(e, g);

		if(type(v_to, g) == NODE && !is_in(v_to, vertices)) {
			auto& m = *matrix(v_to, g);
			if(m.type == F || m.type == DIAGONAL) {

				auto mergeable_vertices = path_closing(gbn, std::vector<Vertex>({v, v_to}));
				for(auto vertex : mergeable_vertices)
					if(type(vertex, g) ==  NODE) {
						auto& m_vertex = *matrix(vertex,g);

						if(!(m_vertex.type == F || m_vertex.type == DIAGONAL)) path_closed = false;
					}

				if(path_closed) {
					auto new_vertices = find_recursively_out(gbn, v_to);
					vertices.insert(new_vertices.begin(), new_vertices.end());
				}
			}
		}
	}

	return vertices;
}

std::set<Vertex> find_recursively_in(GBN gbn, Vertex v) { //can there be loops in a gbn?
    auto g = gbn.graph;
    bool path_closed = true;

    std::set<Vertex> vertices;
    vertices.insert(v);

    for(auto e : boost::make_iterator_range(boost::in_edges(v, g))) {

        Vertex v_from = boost::source(e, g);

        if(type(v_from, g) == NODE && !is_in(v_from, vertices)) {
            auto& m = *matrix(v_from, g);
            if(m.type == F || m.type == DIAGONAL) {

                auto mergeable_vertices = path_closing(gbn, std::vector<Vertex>({v, v_from}));
                for(auto vertex : mergeable_vertices)
                    if(type(vertex, g) ==  NODE) {
                        auto& m_vertex = *matrix(vertex,g);

                        if(!(m_vertex.type == F || m_vertex.type == DIAGONAL)) path_closed = false;
                    }

                if(path_closed) {
                    auto new_vertices = find_recursively_out(gbn, v_from);
                    vertices.insert(new_vertices.begin(), new_vertices.end());
                }
            }
        }
    }

    return vertices;
}

bool reduce_diagonal_matrix(GBN& gbn, Vertex v, std::string& op) {

	op = "reduce_diagonal_matrix";

	auto& g = gbn.graph;
	auto& m = *matrix(v,g);
	std::vector<std::size_t> independent_vars;

	if(m.type != DIAGONAL)
		return false;

	auto& m_DIAG = dynamic_cast<DiagonalMatrix&>(m);

    if(m_DIAG.data.size() % 2 != 0)
        return false;

	//if Diagonal Matrix == ZeroMatrix TODO??: I could also just REPLACE Diagonal for Zero Matrix
	if(m_DIAG.data.empty()) {

		for(auto e_in : boost::make_iterator_range(boost::in_edges(v, g))) {

			Vertex source_vec = boost::source(e_in,g);
			std::size_t source_port = port_from(e_in,g);

			auto v_T = add_vertex(gbn, std::make_shared<TerminatorMatrix>(), "T");
			auto v_Zero = add_vertex(gbn, std::make_shared<ZeroMatrix>(0,1), "zero");

			auto edge_T = boost::add_edge(source_vec, v_T, g).first;
			put(edge_position, g, edge_T, std::make_pair(source_port, 0));
            put(edge_equivalence_class, g, edge_T, equivalence_class(e_in, g));

			for(auto e_out : boost::make_iterator_range(boost::out_edges(v, g))) {
				if(port_from(e_out,g) == port_to(e_in,g)) {
					Vertex target_vec = boost::target(e_out, g);
					std::size_t target_port = port_to(e_out, g);

					auto edge_Zero = boost::add_edge(v_Zero, target_vec, g).first;
					put(edge_position, g, edge_Zero, std::make_pair(0, target_port));
                    put(edge_equivalence_class, g, edge_Zero, equivalence_class(e_out, g));
				}
			}
		}

		remove_vertex(v, gbn);

		return true;
	}

	for(std::size_t var_getting_checked = 0; var_getting_checked < m.n; var_getting_checked++) {
		BitVec assignment(0);
		bool independent = true;

		auto max_i_var = m.n;
		std::size_t i_var = 0;

		while(i_var < max_i_var)  {

			assignment[i_var].flip();
			if(!assignment.test(i_var) && !(assignment.none() && i_var == (max_i_var-1))) {
				i_var++;
			} else {
				auto res1 = m.get(assignment,assignment);
				assignment[var_getting_checked].flip();
				auto res2 = m.get(assignment,assignment);

				if(res1 != res2) {
					independent = false;
					break; //for efficiency
				}
				assignment[var_getting_checked].flip();
				i_var = (assignment.any()) ? 0 : i_var+1;
			}
		}

		if(independent)
			independent_vars.push_back(var_getting_checked);
	}

	if(independent_vars.size() > 0) {
        std::size_t new_matrix_size = m.n - independent_vars.size();

        std::vector<std::size_t> independent_ports = independent_vars;
        /*for(auto independent_var : independent_vars) {
            independent_ports.push_back(m.n - independent_var - 1);
        }*/

        //connect independent ports
        for(auto independent_port : independent_ports) {

            for(auto e_in : boost::make_iterator_range(boost::in_edges(v, g))) {
                if(port_to(e_in,g) == independent_port) {
                    Vertex source_vec = boost::source(e_in,g);
                    std::size_t source_port = port_from(e_in,g);

                    for(auto e_out : boost::make_iterator_range(boost::out_edges(v, g))) {
                        if(port_from(e_out,g) == independent_port) {
                            Vertex target_vec = boost::target(e_out, g);
                            std::size_t target_port = port_to(e_out, g);

                            auto edge_new = boost::add_edge(source_vec, target_vec, g).first;
                            put(edge_position, g, edge_new, std::make_pair(source_port, target_port));
                            put(edge_equivalence_class, g, edge_new, equivalence_class(e_in, g));
                        }
                    }
                }
            }
        }

		//create vertex w/o independent vars
        if(new_matrix_size > 0) {
            auto m_new = std::make_shared<DiagonalMatrix>(new_matrix_size);

            BitVec assignment_old(0);
            BitVec assignment_new(0);
            std::size_t i_var = 0;
            while(i_var < m.n) {
				assignment_old[i_var].flip();
				if (!assignment_old.test(i_var) && !(assignment_old.none() && i_var == (m.n-1))) {
					i_var++;
				} else {
					auto val = m.get(assignment_old, assignment_old);

					std::size_t i_new = 0;
					for(std::size_t i = 0; i < m.n; i++) {
						if(!is_in(i, independent_vars)) {
							assignment_new[i_new] = assignment_old[i];
							i_new++;
						}
					}

					m_new->set(assignment_new, assignment_new, val);
					assignment_new = 0;
					i_var = (assignment_old.any()) ? 0 : i_var+1;
				}
			}

            auto v_new = add_vertex(gbn, m_new, "diagonal");

            for(auto e_in : boost::make_iterator_range(boost::in_edges(v, g))) {
                if(!is_in(port_to(e_in, g), independent_ports)) {
                    Vertex source_vec = boost::source(e_in,g);
                    std::size_t source_port = port_from(e_in,g);
					std::size_t port_diff = std::count_if(independent_ports.begin(), independent_ports.end(), [e_in, g](std::size_t i){return i < port_to(e_in, g);});
					std::size_t target_port = port_to(e_in, g) - port_diff;

					auto edge_new = boost::add_edge(source_vec, v_new, g).first;
                    put(edge_position, g, edge_new, std::make_pair(source_port, target_port));
                    put(edge_equivalence_class, g, edge_new, equivalence_class(e_in, g));
                }
            }

            for(auto e_out : boost::make_iterator_range(boost::out_edges(v, g))) {
                if(!is_in(port_from(e_out, g), independent_ports)) {
                    Vertex target_vec = boost::target(e_out,g);
                    std::size_t target_port = port_to(e_out,g);
                    std::size_t port_diff = std::count_if(independent_ports.begin(), independent_ports.end(), [e_out, g](std::size_t i){return i < port_from(e_out, g);});
                    std::size_t source_port = port_from(e_out, g) - port_diff;

                    auto edge_new = boost::add_edge(v_new, target_vec, g).first;
                    put(edge_position, g, edge_new, std::make_pair(source_port, target_port));
                    put(edge_equivalence_class, g, edge_new, equivalence_class(e_out, g));
                }
            }
        }

        remove_vertex(v, gbn);
        return true;

    } else {
		return false;
	}
}