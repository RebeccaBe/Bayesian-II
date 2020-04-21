#pragma once

#include "../general/gbn.h"
#include "wire_structure.h"
#include "../general/path_closing.h"
#include "../modification/merging.h"



std::vector<Vertex> flip_wire(Wire &wire);

MatrixPtr evaluate(const GBN &gbn);
MatrixPtr evaluate_stepwise(const GBN &gbn);


void update_dependent_wires(WireStructure &wire_structure, Wire independent_wire, std::vector<Vertex> &affected_vertices_vec);
GBN node_elimination (GBN& gbn, std::vector<Vertex> nodes);


