#include "operations_on_gbn.h"
#include "../gbn/modification/vertex_add_remove.h"
#include "../gbn/general/check.h"

#include <iostream>
#include <fstream>

void set_op(const std::vector<std::size_t> places, bool b, GBN& gbn)
{
	auto& g = gbn.graph;
	auto& output_vertices = ::output_vertices(gbn);

	for(auto p : places)
	{
		// get precessor of output port
		auto tmp_in_edges = boost::in_edges(output_vertices[p], g);
		if(std::distance(tmp_in_edges.first, tmp_in_edges.second) != 1)
			throw std::logic_error(std::string("Place ") + std::to_string(p) + " has none or more than one precessor");
		auto e_pre = *(tmp_in_edges.first);
		auto e_pos = get(edge_position, g, e_pre);
		auto v_pre = boost::source(e_pre, g);

		auto v_term = add_vertex(gbn, std::make_shared<TerminatorMatrix>(), "T");
		auto e_term = boost::add_edge(v_pre, v_term, g).first;
		put(edge_position, g, e_term, std::pair<std::size_t, std::size_t>{ e_pos.first, 0 });

		auto v_one = add_vertex(gbn, std::make_shared<OneBMatrix>(b), std::string("1_")+std::to_string(b));
		auto e_one = boost::add_edge(v_one, output_vertices[p], g).first;
		put(edge_position, g, e_one, std::pair<std::size_t, std::size_t>{ 0, 0 });

		boost::remove_edge(v_pre, output_vertices[p], g);
	}
}

void assert_op(const std::vector<std::size_t> places, bool b, GBN& gbn)
{
	auto& g = gbn.graph;
	auto& output_vertices = ::output_vertices(gbn);

	for(auto p : places)
	{
		// get precessor of output port
		auto tmp_in_edges = boost::in_edges(output_vertices[p], g);
		if(std::distance(tmp_in_edges.first, tmp_in_edges.second) != 1)
			throw std::logic_error(std::string("Place ") + std::to_string(p) + " has none or more than one precessor");

		auto e_pre = *(tmp_in_edges.first);
		auto e_pos = get(edge_position, g, e_pre);
		auto v_pre = boost::source(e_pre, g);

		auto v_assert = add_vertex(gbn, std::make_shared<FMatrix>(1,1-b), std::string("F_{1,") + std::to_string(1-b) + "}");
		auto e_term = boost::add_edge(v_pre, v_assert, g).first;
		put(edge_position, g, e_term, std::pair<std::size_t, std::size_t>{ e_pos.first, 0 });
		auto e_out = boost::add_edge(v_assert, output_vertices[p], g).first;
		put(edge_position, g, e_out, std::pair<std::size_t, std::size_t>{ 0, 0 });
		boost::remove_edge(v_pre, output_vertices[p], g);
	}
}

void nassert_op(const std::vector<std::size_t> places, bool b, GBN& gbn)
{
	auto& g = gbn.graph;
	auto& output_vertices = ::output_vertices(gbn);

	std::size_t k = places.size();
	
	auto v_nassert = add_vertex(gbn, std::make_shared<FMatrix>(k,b),std::string("F_{") + std::to_string(k) + "," + std::to_string(b) + "}");
	std::size_t i_place = 0;
	for(auto p : places)
	{
		// get precessor of output port
		auto tmp_in_edges = boost::in_edges(output_vertices[p], g);
		if(std::distance(tmp_in_edges.first, tmp_in_edges.second) != 1)
			throw std::logic_error(std::string("Place ") + std::to_string(p) + " has none or more than one precessor");

		auto e_pre = *(tmp_in_edges.first);
		auto e_pos = get(edge_position, g, e_pre);
		auto v_pre = boost::source(e_pre, g);

		auto e_to_nassert = boost::add_edge(v_pre, v_nassert, g).first;
		put(edge_position, g, e_to_nassert, std::pair<std::size_t, std::size_t>{ e_pos.first, i_place });
		auto e_from_nassert = boost::add_edge(v_nassert, output_vertices[p], g).first;
		put(edge_position, g, e_from_nassert, std::pair<std::size_t, std::size_t>{ i_place, 0 });
		boost::remove_edge(v_pre, output_vertices[p], g);

		i_place++;
	}
}

void setp_op(const std::vector<std::size_t> places, double p, GBN& gbn) //p represents value for ->1
{
	auto& g = gbn.graph;
	auto& output_vertices = ::output_vertices(gbn);

	for(auto p : places)
	{
		// get precessor of output port
		auto tmp_in_edges = boost::in_edges(output_vertices[p], g);
		if(std::distance(tmp_in_edges.first, tmp_in_edges.second) != 1)
			throw std::logic_error(std::string("Place ") + std::to_string(p) + " has none or more than one precessor");
		auto e_pre = *(tmp_in_edges.first);
		auto e_pos = get(edge_position, g, e_pre);
		auto v_pre = boost::source(e_pre, g);

		auto v_term = add_vertex(gbn, std::make_shared<TerminatorMatrix>(), "T");
		auto e_term = boost::add_edge(v_pre, v_term, g).first;
		put(edge_position, g, e_term, std::pair<std::size_t, std::size_t>{ e_pos.first, 0 });

		auto matrix = std::make_shared<DynamicMatrix>(0,1);
		matrix->set(1,0,p);
		matrix->set(0,0,(1-p));
		auto v_p = add_vertex(gbn, matrix, std::to_string(p)+std::string(",")+std::to_string((1-p)));
		auto e_p = boost::add_edge(v_p, output_vertices[p], g).first;
		put(edge_position, g, e_p, std::pair<std::size_t, std::size_t>{ 0, 0 });

		boost::remove_edge(v_pre, output_vertices[p], g);
	}
}

void successp_op(const std::vector<std::vector<std::size_t>> pre_places, const std::vector<std::vector<std::size_t>> post_places,
                 const std::vector<double> probabilities, GBN& gbn) {
    double sum = 0;
    for(auto p : probabilities) {
        sum += p;}
    if(sum-1 > 0.001 || 1-sum > 0.001)
        throw std::logic_error(std::string("Probabilities do not add up to 1. Instead the sum equals " + std::to_string(sum)));

    std::size_t incorporated_transitions = 0;
    std::vector<std::vector<std::size_t>> new_pre_places;
    std::vector<std::vector<std::size_t>> new_post_places;
    std::vector<double> new_probabilities;

    //simplification measure
    for(std::size_t i_transition = 0; i_transition < probabilities.size(); i_transition++) {
		bool valid_pre = validate_transition(pre_places[i_transition], true, gbn);
        if(valid_pre) {
            new_pre_places.push_back(pre_places[i_transition]);
			new_post_places.push_back(post_places[i_transition]);
			new_probabilities.push_back(probabilities[i_transition]);
            incorporated_transitions++;
        }
    }

    if(incorporated_transitions == 0)
        return;
    else if(incorporated_transitions == 1) {
        assert_op(new_pre_places[0],1, gbn);
        check_gbn_integrity(gbn);

        set_op(new_pre_places[0],0,gbn);
        check_gbn_integrity(gbn);

        set_op(new_post_places[0],1,gbn);
        check_gbn_integrity(gbn);
        return;
    }

    std::set<std::size_t> all_places;
	for(auto pre : new_pre_places)
	    all_places.insert(pre.begin(), pre.end());
	for(auto post : new_post_places)
	    all_places.insert(post.begin(), post.end());

	std::map<std::size_t, std::size_t> mapping_place_key;
	std::size_t index = 0;
	for(auto place : all_places) {
	    mapping_place_key.insert(std::pair<std::size_t, std::size_t>(place, index));
	    index++;
	}

	double normalize = 0;
	for(auto p : new_probabilities)
	    normalize += p;

    std::size_t size_n = all_places.size();
	auto matrix = std::make_shared<DynamicMatrix>(size_n, size_n);
	matrix->is_stochastic = false;

	for(std::size_t i = 0; i < new_probabilities.size(); i++){
	    double p = new_probabilities[i]/normalize;
        BitVec assignmentX(0),assignmentY(0);
        std::vector<std::size_t> involved_wires_pre, involved_wires_post;

        for(auto pre : new_pre_places[i]) {
            assignmentY.flip(mapping_place_key.at(pre)); //(X |--pre=1 -- Y=0)
            involved_wires_pre.push_back(mapping_place_key.at(pre));
        }
        for(auto post : new_post_places[i]) {
            assignmentX.flip(mapping_place_key.at(post)); //(X=0 -- post=1 --| Y)
            involved_wires_post.push_back(mapping_place_key.at(post));
        }

        matrix->add(assignmentX,assignmentY,p);

        std::size_t y_wire = 0;
        while (y_wire < size_n){
            if(is_in(y_wire, involved_wires_pre)) { //involved wires will get skipped (we don't want them to have any other values than 0!)
                y_wire++;
                continue;
            } else {
                assignmentY.flip(y_wire);
            }

            if (!assignmentY.test(y_wire)) y_wire++;
            else {
                for(std::size_t x_wire = 0; x_wire < size_n; x_wire++) {
                    if(is_in(x_wire, involved_wires_pre) && !is_in(x_wire, involved_wires_post))
                        assignmentX.set(x_wire, false);
                    else if (is_in(x_wire, involved_wires_post))
                        assignmentX.set(x_wire, true);
                    else
                        assignmentX[x_wire]=assignmentY[x_wire]; //non-involved wires don't change, hence X=Y
                }
                matrix->add(assignmentX,assignmentY,p);
                y_wire = 0;
            }
        }
	}

	auto& g = gbn.graph;
	auto& output_vertices = ::output_vertices(gbn);

	auto v_matrix = add_vertex(gbn, matrix, std::string("successp"));
    std::size_t i_place = 0;
    for(auto p : all_places) {
        // get predecessor of output port
        auto tmp_in_edges = boost::in_edges(output_vertices[p], g);
        if(std::distance(tmp_in_edges.first, tmp_in_edges.second) != 1)
            throw std::logic_error(std::string("Place ") + std::to_string(p) + " has none or more than one precessor");

        auto e_pre = *(tmp_in_edges.first);
        auto e_pos = get(edge_position, g, e_pre);
        auto v_pre = boost::source(e_pre, g);

        auto e_to_matrix = boost::add_edge(v_pre, v_matrix, g).first;
        put(edge_position, g, e_to_matrix, std::pair<std::size_t, std::size_t>{ e_pos.first, i_place });
        auto e_from_matrix = boost::add_edge(v_matrix, output_vertices[p], g).first;
        put(edge_position, g, e_from_matrix, std::pair<std::size_t, std::size_t>{ i_place, 0 });
        boost::remove_edge(v_pre, output_vertices[p], g);

        i_place++;
	}
}



void failp_op(const std::vector<std::vector<std::size_t>> pre_places, const std::vector<double> probabilities, GBN& gbn) {
    double sum = 0;
    for(auto p : probabilities) {
        sum += p;}
    if(sum-1 > 0.001 || 1-sum > 0.001)
        throw std::logic_error(std::string("Probabilities do not add up to 1. Instead the sum equals " + std::to_string(sum)));

    std::set<std::size_t> all_pre_places;
    for(auto pre : pre_places)
        all_pre_places.insert(pre.begin(), pre.end());

    std::map<std::size_t, std::size_t> mapping_place_key;
    std::size_t index = 0;
    for(auto place : all_pre_places) {
        mapping_place_key.insert(std::pair<std::size_t, std::size_t>(place, index));
        index++;
    }

    std::size_t size_n = all_pre_places.size();
    auto matrix = std::make_shared<DiagonalMatrix>(size_n);
    matrix->is_stochastic = false;

    for(std::size_t i = 0; i < probabilities.size(); i++){
        double p = probabilities[i];
        BitVec assignment(0);
        std::vector<std::size_t> involved_wires;

        for(auto pre : pre_places[i])
            involved_wires.push_back(mapping_place_key.at(pre));

        matrix->add(assignment, assignment, p); //(0...0|0...0)

        std::size_t wire = 0;
        while (wire < size_n){
            if(is_in(wire, involved_wires)) { //involved wires will get skipped (we don't want them to have any other values than 0!)
                wire++;
                continue;
            } else {
                assignment.flip(wire);
            }

            if (!assignment.test(wire)) wire++;
            else {
                matrix->add(assignment,assignment,p);
                wire = 0;
            }
        }
    }

    auto& g = gbn.graph;
    auto& output_vertices = ::output_vertices(gbn);

    auto v_matrix = add_vertex(gbn, matrix, std::string("failp"));
    std::size_t i_place = 0;
    for(auto p : all_pre_places) {
        // get predecessor of output port
        auto tmp_in_edges = boost::in_edges(output_vertices[p], g);
        if(std::distance(tmp_in_edges.first, tmp_in_edges.second) != 1)
            throw std::logic_error(std::string("Place ") + std::to_string(p) + " has none or more than one precessor");

        auto e_pre = *(tmp_in_edges.first);
        auto e_pos = get(edge_position, g, e_pre);
        auto v_pre = boost::source(e_pre, g);

        auto e_to_matrix = boost::add_edge(v_pre, v_matrix, g).first;
        put(edge_position, g, e_to_matrix, std::pair<std::size_t, std::size_t>{ e_pos.first, i_place });
        auto e_from_matrix = boost::add_edge(v_matrix, output_vertices[p], g).first;
        put(edge_position, g, e_from_matrix, std::pair<std::size_t, std::size_t>{ i_place, 0 });
        boost::remove_edge(v_pre, output_vertices[p], g);

        i_place++;
    }
}

void successStoch_op(const std::vector<std::vector<std::size_t>> pre_places, const std::vector<std::vector<std::size_t>> post_places,
                 const std::vector<double> probabilities, GBN& gbn) {
    double sum = 0;
    for(auto p : probabilities) {
        sum += p;}
    if(sum-1 > 0.001 || 1-sum > 0.001)
        throw std::logic_error(std::string("Probabilities do not add up to 1. Instead the sum equals " + std::to_string(sum)));

    std::size_t incorporated_transitions = 0;
    std::vector<std::vector<std::size_t>> new_pre_places;
    std::vector<std::vector<std::size_t>> new_post_places;
    std::vector<double> new_probabilities;

    //simplification measure
    for(std::size_t i_transition = 0; i_transition < probabilities.size(); i_transition++) {
        bool valid_pre = validate_transition(pre_places[i_transition], true, gbn);
        if(valid_pre) {
            new_pre_places.push_back(pre_places[i_transition]);
            new_post_places.push_back(post_places[i_transition]);
            new_probabilities.push_back(probabilities[i_transition]);
            incorporated_transitions++;
        }
    }

    if(incorporated_transitions == 0)
        return;
    else if(incorporated_transitions == 1) {
        assert_op(new_pre_places[0],1, gbn);
        check_gbn_integrity(gbn);

        set_op(new_pre_places[0],0,gbn);
        check_gbn_integrity(gbn);

        set_op(new_post_places[0],1,gbn);
        check_gbn_integrity(gbn);
        return;
    }

    std::set<std::size_t> all_places;
    for(auto pre : new_pre_places)
        all_places.insert(pre.begin(), pre.end());
    for(auto post : new_post_places)
        all_places.insert(post.begin(), post.end());

    std::map<std::size_t, std::size_t> mapping_place_key; //maps which place has which wire position
    std::size_t index = 0;
    for(auto place : all_places) {
        mapping_place_key.insert(std::pair<std::size_t, std::size_t>(place, index));
        index++;
    }

    double normalize = 0;
    for(auto p : new_probabilities)
        normalize += p;

    std::size_t size_n = all_places.size();
    auto matrix = std::make_shared<DynamicMatrix>(size_n, size_n);
    matrix->is_stochastic = false;

    for(std::size_t i = 0; i < new_probabilities.size(); i++){
        double p = new_probabilities[i]/normalize;
        BitVec assignmentX(0),assignmentY(0);
        std::vector<std::size_t> involved_wires_pre, involved_wires_post;

        for(auto pre : new_pre_places[i]) {
            assignmentY.flip(mapping_place_key.at(pre)); //(X |--pre=1 -- Y=0)
            involved_wires_pre.push_back(mapping_place_key.at(pre));
        }
        for(auto post : new_post_places[i]) {
            assignmentX.flip(mapping_place_key.at(post)); //(X=0 -- post=1 --| Y)
            involved_wires_post.push_back(mapping_place_key.at(post));
        }

        double probability_norm = p;
        matrix->add(assignmentX,assignmentY,probability_norm);

        std::size_t y_wire = 0;
        while (y_wire < size_n){
            if(is_in(y_wire, involved_wires_pre)) {
                y_wire++;
                continue;
            } else {
                assignmentY.flip(y_wire);
            }

            if (!assignmentY.test(y_wire)) y_wire++;
            else {
                for(std::size_t x_wire = 0; x_wire < size_n; x_wire++) {
                    if(is_in(x_wire, involved_wires_pre) && !is_in(x_wire, involved_wires_post))
                        assignmentX.set(x_wire, false);
                    else if (is_in(x_wire, involved_wires_post))
                        assignmentX.set(x_wire, true);
                    else
                        assignmentX[x_wire]=assignmentY[x_wire]; //non-involved wires don't change, hence X=Y
                }
                double probability_norm = p;
                matrix->add(assignmentX,assignmentY,probability_norm);
                y_wire = 0;
            }
        }
    }

    normalize_matrix_cols(matrix);

    auto& g = gbn.graph;
    auto& output_vertices = ::output_vertices(gbn);

    auto v_matrix = add_vertex(gbn, matrix, std::string("successp"));
    std::size_t i_place = 0;
    for(auto p : all_places) {
        // get predecessor of output port
        auto tmp_in_edges = boost::in_edges(output_vertices[p], g);
        if(std::distance(tmp_in_edges.first, tmp_in_edges.second) != 1)
            throw std::logic_error(std::string("Place ") + std::to_string(p) + " has none or more than one precessor");

        auto e_pre = *(tmp_in_edges.first);
        auto e_pos = get(edge_position, g, e_pre);
        auto v_pre = boost::source(e_pre, g);

        auto e_to_matrix = boost::add_edge(v_pre, v_matrix, g).first;
        put(edge_position, g, e_to_matrix, std::pair<std::size_t, std::size_t>{ e_pos.first, i_place });
        auto e_from_matrix = boost::add_edge(v_matrix, output_vertices[p], g).first;
        put(edge_position, g, e_from_matrix, std::pair<std::size_t, std::size_t>{ i_place, 0 });
        boost::remove_edge(v_pre, output_vertices[p], g);

        i_place++;
    }
}

void normalize_matrix_cols(MatrixPtr matrix) {
    if(matrix->is_stochastic)
        return;

    unsigned long long i_max_row = 1;
    unsigned long long i_max_col = 1;
    i_max_col = i_max_col << matrix->n;
    i_max_row = i_max_row << matrix->m;

    for(unsigned long long i_col = 0; i_col < i_max_col; i_col++) {
        double col_sum = 0;
        for(unsigned long long i_row = 0; i_row < i_max_row; i_row++)
            col_sum += matrix->get(i_row, i_col);

        if (col_sum > 0)
            for(unsigned long long i_row = 0; i_row < i_max_row; i_row++) {
                double old_val = matrix->get(i_row, i_col);
                if(old_val > 0)
                    matrix->set(i_row, i_col, old_val / col_sum);
            }
    }
}

void normalize_matrix_rows(MatrixPtr matrix) {
    unsigned long long i_max_row = 1;
    unsigned long long i_max_col = 1;
    i_max_col = i_max_col << matrix->n;
    i_max_row = i_max_row << matrix->m;

    for(unsigned long long i_row = 0; i_row < i_max_row; i_row++) {
        double row_sum = 0;
        for(unsigned long long i_col = 0; i_col < i_max_col; i_col++)
            row_sum += matrix->get(i_row, i_col);

        if (row_sum > 0)
            for(unsigned long long i_col = 0; i_col < i_max_col; i_col++) {
                double old_val = matrix->get(i_row, i_col);
                if(old_val > 0)
                    matrix->set(i_row, i_col, old_val / row_sum);
            }
    }
}

bool validate_transition(std::vector<std::size_t> places, bool condition, GBN& gbn) {
	auto& g = gbn.graph;
	auto& output_vertices = ::output_vertices(gbn);

	for(auto p : places) {
		auto tmp_in_edges = boost::in_edges(output_vertices[p], g);
		if(std::distance(tmp_in_edges.first, tmp_in_edges.second) != 1)
			throw std::logic_error(std::string("Place ") + std::to_string(p) + " has none or more than one precessor");
		auto v_pre = boost::source(*(tmp_in_edges.first), g);
		auto matrix_v_pre = matrix(v_pre, g);
		if(matrix_v_pre->type == ONE_B) {
			auto& one_matrix = dynamic_cast<OneBMatrix&>(*matrix(v_pre,g));
			if (one_matrix.b != condition){
				return false;
			}
		}
	}
	return true;
}