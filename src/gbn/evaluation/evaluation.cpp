#include "evaluation.h"
#include "../../helpers.hpp"
#include <iostream>
#include "probability_bookkeeper.h"


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
        m = std::make_shared<DiagonalMatrix>(gbn.n, std::vector<BitVec>());
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
    //std::cout << "number of independent wires:" << std::count_if(wires.begin(), wires.end(), [](Wire w){return w.independent;}) << std::endl;
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

MatrixPtr evaluate_stepwise(const GBN &gbn) {
    auto gbn_result = gbn;

    bool changed = false;
    do {
        std::vector<Vertex> nodes;
        auto &g = gbn_result.graph;

        for (auto v : inside_vertices(gbn_result)) {

            bool has_successors = false;
            for (auto e : boost::make_iterator_range(boost::out_edges(v, g))) {
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
            changed = true;
        } else {
            changed = false;
            break;
        }
    }
    while (changed);

    auto vertices = inside_vertices(gbn_result);
    if(vertices.size() > 1) {
        return evaluate(gbn_result);
    } else {
        return matrix(vertices.at(0), gbn_result.graph);
    }
}

GBN node_elimination (GBN& gbn, std::vector<Vertex> nodes) {
    auto &g = gbn.graph;

    std::size_t min_degree = degree(nodes.at(0), g);
    std::size_t max_degree = 0;
    Vertex min_degree_vertex = nodes.at(0);
    for(auto v : nodes) {
        auto v_degree = degree(v, g);
        if(v_degree < min_degree) {
            min_degree = v_degree;
            min_degree_vertex = v;
        }
        if(v_degree > max_degree) {
            max_degree = v_degree;
        }
    }

    std::size_t min_degree_neighbor =  max_degree+1;
    Vertex min_degree_neighbor_vertex = min_degree_vertex;
    for(auto e : boost::make_iterator_range(boost::out_edges(min_degree_vertex, g))) {
        auto successor = boost::target(e, g);
        if(type(successor, g) == NODE) {
            std::vector<Vertex> neighbors = {min_degree_vertex, successor};
            auto path = path_closing(gbn, neighbors);

            if(degree(successor,g) < min_degree_neighbor && path.size()==2) {
                min_degree_neighbor = degree(successor,g);
                min_degree_neighbor_vertex = successor;
            }
        }
    }
    for(auto e : boost::make_iterator_range(boost::in_edges(min_degree_vertex, g))) {
        auto predecessor = boost::source(e, g);
        if(type(predecessor, g) == NODE) {
            std::vector<Vertex> neighbors = {min_degree_vertex, predecessor};
            auto path = path_closing(gbn, neighbors);

            if (degree(predecessor, g) < min_degree_neighbor && path.size() == 2) {
                min_degree_neighbor = degree(predecessor, g);
                min_degree_neighbor_vertex = predecessor;
            }
        }
    }

    std::vector<Vertex> minimal_degree_neighbors = {min_degree_vertex, min_degree_neighbor_vertex};
    merge_vertices(gbn, minimal_degree_neighbors);

    return gbn;
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
