#include "gbn.h"

GBN::GBN(Index n, Index m, Index n_inside_vertices)
	:n(n), m(m), n_vertices(n_inside_vertices + n + m), graph(n_vertices), n_initial_inside_vertices(n_inside_vertices), n_initial_n(n), n_initial_m(m),max_eq_class_counter(0)
{
	for(Index i = 0; i < n_vertices; i++)
	{
		put(vertex_type, graph, i, NODE);
		visible_vertices.insert(i);
	}

	for(Index v = 0; v < n_inside_vertices; v++)
		put(vertex_name, graph, v, std::string("v_")+std::to_string(v));

	for(Index i = 0; i < n; i++)
	{
		put(vertex_type, graph, n_inside_vertices+i, INPUT);
		put(vertex_name, graph, n_inside_vertices+i, std::string("i_")+std::to_string(i));
		input_vertices.push_back(n_inside_vertices+i);
	}
	for(Index i = 0; i < m; i++)
	{
		put(vertex_type, graph, n_inside_vertices+n+i, OUTPUT);
		put(vertex_name, graph, n_inside_vertices+n+i, std::string("o_")+std::to_string(i));
		output_vertices.push_back(n_inside_vertices+n+i);
	}
}

Index n(const GBNGraph::vertex_descriptor& v, const GBNGraph& g)
{
	return get(vertex_matrix, g, v)->n;
}

Index m(const GBNGraph::vertex_descriptor& v, const GBNGraph& g)
{
	return get(vertex_matrix, g, v)->m;
}


Index port_from(const GBNGraph::edge_descriptor& e, const GBNGraph& g)
{
	return get(edge_position, g, e).first;
}

Index port_to(const GBNGraph::edge_descriptor& e, const GBNGraph& g)
{
	return get(edge_position, g, e).second;
}

Index input_idx(const Vertex& v, const GBN& gbn)
{
	return v-gbn.n_initial_inside_vertices;
}

Index output_idx(const Vertex& v, const GBN& gbn)
{
	return v-gbn.n_initial_inside_vertices-gbn.n_initial_n;
}

const std::string& name(const GBNGraph::vertex_descriptor& v, const GBNGraph& g)
{
	return get(vertex_name, g, v);
}

VertexType type(const GBNGraph::vertex_descriptor& v, const GBNGraph& g)
{
	return get(vertex_type, g, v);
}

MatrixPtr matrix(const GBNGraph::vertex_descriptor& v, const GBNGraph& g)
{
	return get(vertex_matrix, g, v);
}




std::vector<Vertex> all_vertices(const GBN& gbn)
{
	return std::vector<Vertex>(gbn.visible_vertices.begin(), gbn.visible_vertices.end());
}

std::vector<Vertex> inside_vertices(const GBN& gbn)
{
	std::vector<Vertex> rtn;
	std::copy_if(gbn.visible_vertices.begin(), gbn.visible_vertices.end(), std::back_inserter(rtn), [&gbn](const Vertex& v) {
		return type(v,gbn.graph) == NODE;
	});

	return rtn;
}

const std::vector<Vertex>& input_vertices(const GBN& gbn)
{
	return gbn.input_vertices;
}

const std::vector<Vertex>& output_vertices(const GBN& gbn)
{
	return gbn.output_vertices;
}

std::size_t degree(const GBNGraph::vertex_descriptor& v, const GBNGraph& g) {
    std::size_t degree = 0;
    std::vector<std::size_t> ports;

    degree += boost::make_iterator_range(boost::in_edges(v,g)).size();

    for(auto e : boost::make_iterator_range(boost::out_edges(v,g))) {
        std::size_t port = port_from(e, g);
        if(std::find(ports.cbegin(), ports.cend(), port) == ports.cend()) {
            degree++;
            ports.push_back(port);
        }
    }
    return degree;
}

std::vector<Vertex> neighbors(const GBNGraph::vertex_descriptor& v, const GBNGraph& g) {
    std::vector<Vertex> neighbors;

    for(auto e : boost::make_iterator_range(boost::out_edges(v,g))) {
        auto successor = boost::target(e, g);
        if(type(successor, g) == NODE)
            neighbors.push_back(successor);
    }

    for(auto e : boost::make_iterator_range(boost::in_edges(v,g))) {
        auto predecessor = boost::source(e, g);
        if(type(predecessor, g) == NODE)
            neighbors.push_back(predecessor);
    }

    return neighbors;
}

std::vector<Vertex> all_neighbors(const GBNGraph::vertex_descriptor& v, const GBNGraph& g) {
    std::vector<Vertex> neighbors;

    for(auto e : boost::make_iterator_range(boost::out_edges(v,g))) {
        auto successor = boost::target(e, g);
        neighbors.push_back(successor);
    }

    for(auto e : boost::make_iterator_range(boost::in_edges(v,g))) {
        auto predecessor = boost::source(e, g);
        neighbors.push_back(predecessor);
    }

    return neighbors;
}

std::vector<Edge> all_edges(std::vector<Vertex>& vertices, const GBNGraph& g) {
    std::set<Edge> all_edges;

    for(auto v : vertices) {
        for(auto e : boost::make_iterator_range(boost::out_edges(v,g))) {
                all_edges.insert(e);
        }
        for(auto e : boost::make_iterator_range(boost::in_edges(v,g))) {
                all_edges.insert(e);
        }
    }

    std::vector<Edge> all_edges_vector (all_edges.begin(), all_edges.end());
     return all_edges_vector;

     //ODER:
    //return boost::edges(g);
}

std::size_t equivalence_class (const GBNGraph::edge_descriptor& e, const GBNGraph& g) {
    return get(edge_equivalence_class, g, e);
}

std::vector<Edge> edges_of_equivalence_class(std::size_t eq_class, const GBN& gbn) {
    std::vector<Edge> edges;
    auto vertices = all_vertices(gbn);

    for(auto e : all_edges(vertices, gbn.graph)) {
        if(equivalence_class(e, gbn.graph) == eq_class)
            edges.push_back(e);
    }
    return edges;
}

std::size_t max_equivalence_counter_and_increase(GBN& gbn){
    std::size_t max_equivalence_class = gbn.max_eq_class_counter;
    gbn.max_eq_class_counter++;
    return max_equivalence_class;
}

Edge get_edge(const Port& v_from, const Port& v_to, const GBNGraph& g) {
    for(auto e : boost::make_iterator_range(boost::out_edges(v_from.first,g))) {
        if(boost::target(e, g) == v_to.first && port_from(e,g)==v_from.second && port_to(e,g)==v_to.second)
            return e;
    }
    throw std::logic_error(std::string("Edge does not exist"));
}

Edge get_edge(const Vertex& v_from, const Vertex& v_to, const GBNGraph& g) {
    auto e = boost::edge(v_from, v_to, g);

    if(!e.second)
        throw std::logic_error(std::string("Edge does not exist"));

    return e.first;
}