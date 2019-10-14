#pragma once

#include "../general/gbn.h"
#include "wire_structure.h"


namespace improved {

    std::vector<Vertex> flip_wire(Wire &wire);

    MatrixPtr evaluate(const GBN &gbn);

    void update_dependent_wires(WireStructure &wire_structure, Wire independent_wire, std::vector<Vertex> &affected_vertices_vec);

}