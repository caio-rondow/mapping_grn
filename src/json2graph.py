import json
import networkx as nx
import pydot

f = open('mesh_8x8.json')

data = json.load(f)


# read json and get id and neighbors 
def get_nodes_and_edges(data):
    nodes = []
    edges = {}

    for i in data['PE']:
        node = i['id']
        nodes.append(node)
        edges[node] = i['neighbors']

    return nodes, edges

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


G = make_digraph(data)
PG = nx.nx_pydot.to_pydot(G)

print(PG)
