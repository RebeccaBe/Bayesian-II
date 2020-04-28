#pragma once

#include "../general/gbn.h"
#include "wire_structure.h"
#include "../general/path_closing.h"
#include "../modification/merging.h"
#include "../../cnu/operations_on_gbn.h"



std::vector<Vertex> flip_wire(Wire &wire);

MatrixPtr evaluate(const GBN &gbn);
MatrixPtr evaluate_stepwise(const GBN &gbn);
MatrixPtr evaluate_specific_place(std::size_t place, const GBN &gbn);

