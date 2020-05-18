#include "wire_structure.h"

#include "../../helpers.hpp"
#include <iostream>



WireStructure build_wire_structure(const GBN &gbn) {
    auto &g = gbn.graph;

    auto input_vertices = ::input_vertices(gbn);
    auto output_vertices = ::output_vertices(gbn);
    auto inside_vertices = ::inside_vertices(gbn);
    std::vector<Vertex> origin_vertices;
    std::copy(input_vertices.begin(), input_vertices.end(), std::back_inserter(origin_vertices));
    std::copy(inside_vertices.begin(), inside_vertices.end(), std::back_inserter(origin_vertices));

    // reserve memory for bitvecs
    WireStructure wire_structure;
    auto all_v = all_vertices(gbn);
    auto v_max = *std::max_element(all_v.begin(), all_v.end()); // TODO: put v_max into GBN as a property
    wire_structure.vertex_input_bitvecs.resize(v_max + 1);
    wire_structure.vertex_output_bitvecs.resize(v_max + 1);
    for (auto v : ::inside_vertices(gbn)) {
        wire_structure.vertex_input_bitvecs[v] = std::make_shared<BitVec>();
        wire_structure.vertex_output_bitvecs[v] = std::make_shared<BitVec>();
    }
    wire_structure.input_bitvec = std::make_shared<BitVec>();
    wire_structure.output_bitvec = std::make_shared<BitVec>();

    // build a mapping from all vertex outputs to vertex inputs
    std::map<Port, std::set<Port>> from_to_port_map;
    for (auto v : origin_vertices) {
        for (auto e : boost::make_iterator_range(boost::out_edges(v, g))) {
            Vertex v_from = boost::source(e, g), v_to = boost::target(e, g);
            Port p_from{v_from, port_from(e, g)}, p_to{v_to, port_to(e, g)};
            from_to_port_map[p_from].insert(p_to);
        }
    }

    std::map<Port, std::size_t> master_wires_map;
    std::size_t name = 0; //v1

    // build actual wire structure
    for (auto[p_from, ports_to] : from_to_port_map) {
        Wire w;

        w.name = name;
        name++; //v1
        w.source = p_from;

        if (is_in(p_from.first, input_vertices))
            w.io_ports.push_back({wire_structure.input_bitvec, input_idx(p_from.first, gbn)});
        else
            w.inside_ports.push_back(
                    {p_from.first, wire_structure.vertex_output_bitvecs[p_from.first], p_from.second});

        for (auto p_to : ports_to) {
            if (is_in(p_to.first, output_vertices))
                w.io_ports.push_back({wire_structure.output_bitvec, output_idx(p_to.first, gbn)});
            else
                w.inside_ports.push_back(
                        {p_to.first, wire_structure.vertex_input_bitvecs[p_to.first], p_to.second});
        }

        //build dependency
        if (type(p_from.first, g) == INPUT) {
            master_wires_map.insert(std::pair<Port, std::size_t>(p_from, w.name));

        } else {
            const auto &matrix_p_from = *matrix(p_from.first, g);
            if (!(matrix_p_from.type == F || matrix_p_from.type == DIAGONAL) && matrix_p_from.diag_places.empty())
                master_wires_map.insert(std::pair<Port, std::size_t>(p_from, w.name));
            else if (!matrix_p_from.diag_places.empty()){
                bool is_independent = true;
                for(auto [input, output] : matrix_p_from.diag_places)
                    if(output == p_from.second)
                        is_independent = false;
                if(is_independent)
                    master_wires_map.insert(std::pair<Port, std::size_t>(p_from, w.name));
            }
            else
                w.independent = false;
        }

        wire_structure.wires.push_back(w);
    }


    //build actual dependency
    for (auto &wire : wire_structure.wires) {
        if (!wire.independent) {

            bool still_dependent = true;
            Port deep_source_Port = wire.source;

            while (still_dependent) {

                auto deep_source_Vertex = deep_source_Port.first;
                auto edge_list = boost::make_iterator_range(boost::in_edges(deep_source_Vertex, g));
                for (auto edge : edge_list) {

                    const auto &deep_source_Matrix = *matrix(deep_source_Vertex, g);
                    switch (deep_source_Matrix.type) {
                        case F:
                        case DIAGONAL:
                            if (port_to(edge, g) == deep_source_Port.second) {
                                deep_source_Port = std::pair<Vertex, std::size_t>(boost::source(edge, g),
                                                                                  port_from(edge, g));

                                auto search = master_wires_map.find(deep_source_Port);
                                if (search != master_wires_map.end()) {
                                    wire.master_wire = master_wires_map.at(deep_source_Port);
                                    still_dependent = false;
                                    //std::cout << "wire: " << wire.name << " dependent of: " << wire.master_wire << std::endl;
                                    break;
                                } else break;
                            }
                            break;
                        case DYNAMIC:
                            if (port_to(edge, g) == deep_source_Matrix.diag_places.at(deep_source_Port.second)) {
                                deep_source_Port = std::pair<Vertex, std::size_t>(boost::source(edge, g),
                                                                                  port_from(edge, g));

                                auto search = master_wires_map.find(deep_source_Port);
                                if (search != master_wires_map.end()) {
                                    wire.master_wire = master_wires_map.at(deep_source_Port);
                                    still_dependent = false;
                                    //std::cout << "wire: " << wire.name << " dependent of: " << wire.master_wire << std::endl;
                                    break;
                                } else break;
                            }
                            break;
                    }
                }
            }
        }
    }

    return wire_structure;
}

void print_wire_structure(std::ostream &ostr, const WireStructure &wire_structure, const GBN &gbn) {
    ostr << "n_wires: " << wire_structure.wires.size() << std::endl;

    std::size_t i_wire = 0;
    for (auto w : wire_structure.wires) {
        ostr << "--- " << w.name << " --- " << std::endl;
        ostr << "inside_ports: " << std::endl;
        for (auto[v, bitvec, pos] : w.inside_ports) {
            ostr << name(v, gbn.graph) << " ";
            ostr << pos << " ";
            ostr << std::endl;
        }
        ostr << "io_ports: ";
        for (auto[bitvec, pos] : w.io_ports)
            ostr << pos << " ";
        ostr << std::endl;
        if (!w.independent)
            ostr << "master wire: " << w.master_wire << std::endl;
    }
}

