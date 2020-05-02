#include "evaluation.h"
#include "../../helpers.hpp"
#include <iostream>
#include <fstream>
#include "probability_bookkeeper.h"
#include "../general/gbn_io.h"
#include "../modification/vertex_add_remove.h"
#include "../general/subgbn.h"
#include "../simplification/global_simplification.h"
#include "../simplification/local_simplification.h"


// returns which vertices have to be updated
// if wire. active == false afterwards hints towards carry bit
std::vector<Vertex> flip_wire(Wire &wire) {
    std::vector<Vertex> changed_vertices;
    wire.active = !wire.active;

    for (auto[v, p_bitvec, i_bit] : wire.inside_ports) {
        // TODO: following if could be replaced by one flip without branching?
        if (wire.active)
            p_bitvec->set(i_bit);
        else
            p_bitvec->reset(i_bit);
        changed_vertices.push_back(v); // TODO: could be optimized out by using the vertices from wire...
    }

    for (auto[p_bitvec, i_bit] : wire.io_ports) {
        p_bitvec->flip(i_bit);
        // TODO: following if could be replaced by one flip without branching?
        if (wire.active)
            p_bitvec->set(i_bit);
        else
            p_bitvec->reset(i_bit);
    }

    return changed_vertices;
}

void update_dependent_wires(WireStructure &wire_structure, Wire independent_wire, std::vector<Vertex> &affected_vertices_vec) {
    for (Wire &w : wire_structure.wires) {
        if (!w.independent) {
            if (w.master_wire == independent_wire.name) {
                auto new_affected_vertices_vec = flip_wire(w);
                affected_vertices_vec.insert(affected_vertices_vec.end(), new_affected_vertices_vec.begin(), new_affected_vertices_vec.end());
            }
        }
    }
}

MatrixPtr evaluate(const GBN &gbn) {
    auto &g = gbn.graph;
    auto wire_structure = build_wire_structure(gbn);
    auto &wires = wire_structure.wires;

    bool is_dynamic = 0;
    for(auto vertex : all_vertices(gbn)) {
        if(type(vertex, g) == NODE) {
            auto matrix_type = matrix(vertex, g)->type;
            if (!(matrix_type == F || matrix_type == DIAGONAL))
                is_dynamic = 1;
        }
    }

    MatrixPtr m;
    if(!is_dynamic && (input_vertices(gbn).size() == output_vertices(gbn).size()))
        m = std::make_shared<DiagonalMatrix>(gbn.n);
    else
        m = std::make_shared<DynamicMatrix>(gbn.n, gbn.m);

    auto inside_vertices = ::inside_vertices(gbn);

    // init probabilities to the value at zero
    auto all_v = all_vertices(gbn);
    auto v_max = *std::max_element(all_v.begin(), all_v.end()); // TODO: put v_max into GBN as a property
    ProbabilityBookkeeper bk(v_max + 1, inside_vertices);
    for (auto v : inside_vertices) {
        auto &m_v = *matrix(v, g);
        auto p = m_v.get(*wire_structure.vertex_output_bitvecs[v], *wire_structure.vertex_input_bitvecs[v]);
        bk.update_one_node(v, p);
    }
    double product = bk.get_product();
    m->add(*wire_structure.output_bitvec, *wire_structure.input_bitvec, product);

    std::size_t i_wire = 0;
    std::cout << "number of independent wires:" << std::count_if(wires.begin(), wires.end(), [](Wire w){return w.independent;}) << std::endl;
    const std::size_t max_i_wire = wires.size();

    std::set<Vertex> affected_vertices; // TODO: optimization: replace this with flat set
    while (i_wire < max_i_wire) // TODO: do this more efficiently with gray codes -> only one wire flip at a time needed
    {
        auto &w = wires[i_wire];

        if (w.independent) {

            auto affected_vertices_vec = flip_wire(w);
            update_dependent_wires(wire_structure, w, affected_vertices_vec);
            affected_vertices.insert(affected_vertices_vec.begin(), affected_vertices_vec.end());

            if (!w.active) // carry bit needed
            {
                i_wire++;
            } else {
                for (auto v : affected_vertices) {
                    const auto &m_v = *matrix(v, g);
                    auto p = m_v.get(*wire_structure.vertex_output_bitvecs[v],
                                     *wire_structure.vertex_input_bitvecs[v]);
                    bk.update_one_node(v, p);
                }

                double product = bk.get_product();
                m->add(*wire_structure.output_bitvec, *wire_structure.input_bitvec, product);
                affected_vertices.clear();
                i_wire = 0;
            }
        } else
            i_wire++;
    }

    for (auto v : inside_vertices) {
        auto &m_v = *matrix(v, g);
        if (!m_v.is_stochastic) {
            m->is_stochastic = false;
            break;
        }
    }

    return m;
}

std::size_t get_edges(std::vector<Vertex> vertices, const GBN gbn) {

    auto sub_gbn = SubGBN::make_from_vertices(gbn, vertices);

    auto wire_structure = build_wire_structure(sub_gbn.gbn);
    auto wires = wire_structure.wires;
    return std::count_if(wires.begin(), wires.end(), [](Wire w){return w.independent;});
}

std::vector<Vertex> expand_list(GBN gbn, std::vector<Vertex> nodes) {
    auto g = gbn.graph;

    std::vector<Vertex> expanded_list = nodes;
    std::vector<Vertex> neighbors_v;
    bool changed = false;

    for(auto n : nodes) {
        auto neighbors_n = all_neighbors(n, g);
        neighbors_v.insert(std::end(neighbors_v), std::begin(neighbors_n), std::end(neighbors_n));
    }

    for(auto n : neighbors_v) {
        if(type(n, gbn.graph) == NODE) {

            auto p_m = matrix(n, g);
             if(p_m->type == F || p_m->type == DIAGONAL) {
                 bool all_input_edges_in_set, all_output_edges_in_set = true;

                 for(auto e : boost::make_iterator_range(boost::out_edges(n,g))) {
                     auto successor = boost::target(e, g);
                     if(!is_in(successor, nodes)) {
                         all_input_edges_in_set = false;
                         break;
                     }
                 }

                 if(!all_input_edges_in_set)
                     for(auto e : boost::make_iterator_range(boost::in_edges(n,g))) {
                         auto predecessor = boost::source(e, g);
                         if(!is_in(predecessor, nodes)) {
                             all_output_edges_in_set = false;
                             break;
                         }
                     }

                 if((all_input_edges_in_set || all_output_edges_in_set) && !is_in(n, expanded_list)) {
                     expanded_list.push_back(n);
                     changed = true;
                     continue;
                 }
             }
        }

        auto neighbors_n = all_neighbors(n, g);
        bool neighbors_included_in_nodes = true;

        for(auto nn : neighbors_n) {
            if(!is_in(nn, nodes)) {
                neighbors_included_in_nodes = false;
                break;
            }
        }

        if(neighbors_included_in_nodes && !is_in(n, expanded_list) && type(n, g) == NODE)
            expanded_list.push_back(n);
    }

    if(changed)
        return expand_list(gbn, expanded_list);
    else
        return expanded_list;
}

void find_diagonal_nodes(GBN gbn, std::set<Vertex>& nodes, std::vector<Vertex>& nodes_to_be_explored) {
    auto g = gbn.graph;
    std::vector<Vertex> nodes_to_be_explored_new;

    for(auto n : nodes_to_be_explored) {
        if(type(n, g) == NODE) {
            auto p_m = matrix(n, g);

            if (p_m->type == F || p_m->type == DIAGONAL) nodes.insert(n);
            auto neighbors_n = all_neighbors(n, g);
            for (auto neighbor : neighbors_n) {
                if(type(neighbor, g) == NODE) {
                    auto p_m = matrix(n, g);

                    if (!is_in(neighbor, nodes) && !is_in(neighbor, nodes_to_be_explored_new) &&
                        (p_m->type == F || p_m->type == DIAGONAL)) {
                        nodes_to_be_explored_new.push_back(neighbor);
                    }
                }
            }
        }
    }
    if(nodes_to_be_explored_new.empty()) {
        return;
    } else {
        find_diagonal_nodes(gbn, nodes, nodes_to_be_explored_new);
    }
}

GBN node_elimination(GBN& gbn, std::vector<Vertex> nodes) {
    auto &g = gbn.graph;

    std::size_t min_degree = get_edges(nodes, gbn) + 1;
    std::vector<Vertex> min_degree_neighbors = {nodes[0], nodes[1]}; //dummy values; are these even  necessary?

    bool diagonal_vertices_exist = false;
    //Can I find a case, where I can actually use this for efficiency purposes?
    /*for(auto v : nodes) {
        auto p_m = matrix(v, g);
        if(p_m->type == F || p_m->type == DIAGONAL) {

            std::vector<Vertex> vertex;
            vertex.push_back(v);

            std::set<Vertex> diagonal_vertices;
            find_diagonal_nodes(gbn, diagonal_vertices, vertex);

            if(diagonal_vertices.size() > 1) {
                std::vector<Vertex> diagonal_vertices_vector(diagonal_vertices.begin(), diagonal_vertices.end());
                min_degree_neighbors = diagonal_vertices_vector;
                diagonal_vertices_exist = true;
                for(auto n : min_degree_neighbors) {
                    if(type(n, g) == NODE) {
                        auto p_m = matrix(n, g);
                    }
                }
            }
        }
    }*/

    //First get rid of the Terminator nodes
    bool terminator_exists = false;
    if(!diagonal_vertices_exist) {
        for (auto v : nodes) {
            if (m(v, g) == 0) {
                terminator_exists = true;
                auto neighbors_v = neighbors(v, g); //this checks all vertex pairs at least twice

                for (auto neighbor_v : neighbors_v) {
                    std::vector<Vertex> neighbors_i = {v, neighbor_v};

                    std::size_t v_degree = get_edges(neighbors_i, gbn);
                    if (v_degree < min_degree) {
                        min_degree = v_degree;
                        min_degree_neighbors = neighbors_i;
                    }
                }
            }
        }
    }

    //Second get rid of nodes without input
    bool no_input_exists = false;
    if(!terminator_exists && !diagonal_vertices_exist) {
        for (auto v : nodes) {
            if (n(v, g) == 0) {
                no_input_exists = true;
                auto neighbors_v = neighbors(v, g); //this checks all vertex pairs at least twice

                for (auto neighbor_v : neighbors_v) {
                    std::vector<Vertex> neighbors_i = {v, neighbor_v};

                    std::size_t v_degree = get_edges(neighbors_i, gbn);
                    if (v_degree < min_degree) {
                        min_degree = v_degree;
                        min_degree_neighbors = neighbors_i;
                    }
                }
            }
        }
    }

    //Are their loops between 2 nodes? Eliminate that next
    bool loops_exist = false;
    /*if(!terminator_exists && !no_input_exists && !diagonal_vertices_exist) {
        for(auto v : nodes) {
            if(n(v,g) > 0 && m(v,g) > 0) {
                std::set<Vertex> successors, predecessors;

                for(auto e : boost::make_iterator_range(boost::out_edges(v,g)))
                    successors.insert(boost::target(e, g));

                for(auto e : boost::make_iterator_range(boost::in_edges(v,g)))
                    predecessors.insert(boost::source(e, g));

                for(auto predecessor : predecessors) {
                    std::vector<Vertex> neighbors_i = {v, predecessor};

                    if(is_in(predecessor, successors)) {
                        terminator_exists = true;

                        std::size_t v_degree = get_edges(neighbors_i, gbn);
                        if (v_degree < min_degree) {
                            min_degree = v_degree;
                            min_degree_neighbors = neighbors_i;
                        }
                        break; // for efficiency
                    }
                }
            }
        }
    }*/

    if (!terminator_exists && !no_input_exists && !loops_exist && !diagonal_vertices_exist) {
        for (auto v : nodes) {
            auto neighbors_v = neighbors(v, g); //this checks all vertex pairs at least twice

            for (auto neighbor_v : neighbors_v) {
                std::vector<Vertex> neighbors_i = {v, neighbor_v};

                std::size_t v_degree = get_edges(neighbors_i, gbn);
                if (v_degree < min_degree) {
                    min_degree = v_degree;
                    min_degree_neighbors = neighbors_i;
                }
            }
        }
    }

    auto expanded_nodes = expand_list(gbn, min_degree_neighbors);
    for(auto node : expanded_nodes) {
        if(n(node, g) == 0) {
            std::string op;
            normalize_substoch_front_vertices_without_inputs(gbn, node, op);
        }
    }

    merge_vertices(gbn, expanded_nodes);

    return gbn;
}

void find_relevant_nodes(GBN gbn, std::set<Vertex>& nodes, std::vector<Vertex>& nodes_to_be_explored) {
    auto g = gbn.graph;
    auto nodes_copy = nodes;
    std::vector<Vertex> nodes_to_be_explored_new;

    for(auto n : nodes_to_be_explored) {
        if(type(n, g) == NODE) nodes.insert(n);
        auto neighbors_n  = all_neighbors(n, g);
        for(auto neighbor : neighbors_n) {
            if(!is_in(neighbor, nodes) && !is_in(neighbor, nodes_to_be_explored_new) && type(neighbor, g) == NODE) {
                nodes_to_be_explored_new.push_back(neighbor);
            }
        }
    }
    if(nodes_to_be_explored_new.empty()) {
        return;
    } else {
        find_relevant_nodes(gbn, nodes, nodes_to_be_explored_new);
    }
}

MatrixPtr evaluate_stepwise(const GBN &gbn) {
    auto gbn_result = gbn;
    int index = 0;

    bool changed = false;
    do {
        std::vector<Vertex> nodes;
        auto &g = gbn_result.graph;

        for (auto v : inside_vertices(gbn_result)) {

            bool has_successors = false;
            for (auto e : boost::make_iterator_range(boost::out_edges(v, g))) { //use neighbor function
                if (type(boost::target(e, g), g) != OUTPUT) {
                    has_successors = true;
                    break;
                }
            }

            bool has_predecessors = false;
            for (auto e : boost::make_iterator_range(boost::in_edges(v, g))){
                if (type(boost::source(e, g), g) != INPUT) {
                    has_predecessors = true;
                    break;
                }
            }

            if (has_predecessors || has_successors)
                nodes.push_back(v);
        }

        if(!nodes.empty()) {
            gbn_result  = node_elimination(gbn_result, nodes);

            std::ofstream f1("Evaluation-step"+std::to_string(index)+".dot");
            draw_gbn_graph(f1, gbn_result, "");
            std::ofstream f("Evaluation-step.dot");
            draw_gbn_graph(f, gbn_result, "");
            index++;

            changed = true;
        } else {
            changed = false;
            break;
        }
    }
    while (changed);

    auto vertices = inside_vertices(gbn_result);
    MatrixPtr m;
    if(vertices.size() > 1) {
        m = evaluate(gbn_result);
    } else {
        m = matrix(vertices.at(0), gbn_result.graph);
    }

    normalize_matrix_cols(m);
    return m;
}

MatrixPtr evaluate_specific_place(std::size_t place, const GBN &gbn_og) {
    auto gbn = gbn_og;
    auto &g = gbn.graph;
    auto output_vertices = ::output_vertices(gbn);

    std::vector<Vertex> chosen_output;
    chosen_output.push_back(output_vertices[place]);

    for(std::size_t p = 0; p < gbn.m; p++) {
        if(p == place) continue;

        // get predecessor of output port
        auto tmp_in_edges = boost::in_edges(output_vertices[p], g);
        if(std::distance(tmp_in_edges.first, tmp_in_edges.second) != 1)
            throw std::logic_error(std::string("Place ") + std::to_string(p) + " has none or more than one precessor");
        auto e_pre = *(tmp_in_edges.first);
        auto e_pos = get(edge_position, g, e_pre);
        auto v_pre = boost::source(e_pre, g);

        auto v_term = add_vertex(gbn, std::make_shared<TerminatorMatrix>(), "T");
        auto e_term = boost::add_edge(v_pre, v_term, g).first;
        put(edge_position, g, e_term, std::pair<std::size_t, std::size_t>{ e_pos.first, 0 });

        boost::remove_edge(v_pre, output_vertices[p], g);
    }

    std::set<Vertex> relevant_nodes;
    find_relevant_nodes(gbn, relevant_nodes, chosen_output);
    std::vector<Vertex> relevant_nodes_vector (relevant_nodes.begin(), relevant_nodes.end());

    auto sub_gbn = SubGBN::make_from_vertices(gbn, relevant_nodes_vector);

    std::string s = "";
    //merge_F_matrices_to_diagonal_matrix(sub_gbn.gbn, inside_vertices(sub_gbn.gbn)[30], s);
    std::cout << s << std::endl;

    return evaluate_stepwise(sub_gbn.gbn);
}