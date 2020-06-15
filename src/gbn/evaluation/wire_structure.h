#pragma once

#include "../general/gbn.h"



struct Wire {
    std::vector<std::tuple<Vertex, BitVecPtr, Index>> inside_ports;
    std::vector<std::pair<BitVecPtr, Index>> io_ports;
    bool active = false; // TODO: remove this and only flip bits?
    Port source;

    bool independent = true;
    std::size_t master_wire;
    std::size_t equivalence_class;
    //std::vector<std::size_t> slaveWires; //TODO: anderer Ansatz
    std::size_t name; //v1
};


struct WireStructure {
    std::vector<BitVecPtr> vertex_input_bitvecs;
    std::vector<BitVecPtr> vertex_output_bitvecs;

    BitVecPtr input_bitvec;
    BitVecPtr output_bitvec;

    std::vector<Wire> wires;

    std::map<std::size_t, std::vector<Wire>> equivalence_classes_vec;
};


WireStructure build_wire_structure(const GBN &gbn);

void print_wire_structure(std::ostream &ostr, const WireStructure &wire_structure, const GBN &gbn);