import json
import networkx as nx
import pydot


# read json and get id and neighbors 
def get_nodes_and_edges(data):
    nodes = []
    edges = {}

    for i in data['pe']:
        node = i['id']
        nodes.append(node)
        edges[node] = i['neighbors']

    return nodes, edges

def get_grid_dimension(data):
    return data['data_width']

# build digraph
def make_digraph(data):
    nodes,edges = get_nodes_and_edges(data)

    G = nx.DiGraph()

    for node in nodes:
        G.add_node(node)
    
    for key in edges:
        for i in range(len(edges[key])):
            G.add_edge(key,edges[key][i])

    return G


# networkx graph to .dot
def nx_2_dot(graph):
    return nx.nx_pydot.to_pydot(graph) 
