#pragma once

#include "../general/gbn.h"
#include "wire_structure.h"
#include "../general/path_closing.h"
#include "../modification/merging.h"
#include "../../cnu/operations_on_gbn.h"

enum EvaluationType {
    DEFAULT,
    FILLIN,
    DEGREE
};

std::vector<Vertex> flip_wire(Wire &wire);

MatrixPtr evaluate(const GBN &gbn, EvaluationType eval_type = DEFAULT);
MatrixPtr evaluate_specific_place(std::size_t place, const GBN &gbn, EvaluationType eval_type = DEGREE);

