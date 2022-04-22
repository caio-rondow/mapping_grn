import networkx as nx 
import pydot
import pandas as pd
import matplotlib.pyplot as plt

'''
    1. cgra dimensions
    2. list of nodes
    3. an subgraph with rank=same
'''


def sa_curve(data):
    df = pd.DataFrame(data)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(df[1], df[0], color='tab:blue')
    plt.show()


def print_dot(dot: pydot.Dot) -> None:

    string_dot = dot.to_string()
    string_dot = string_dot.replace("strict digraph  {", "digraph layout  {")


    text_file = open("Output.txt", "w")
    text_file.write(string_dot)
    text_file.close()


def build_dot(graph: pydot.Dot, nodes: list, dim: list) -> None:

    '''
        dim = row x col -> dim[0] = row din[1] = col
    '''

    graph.create_attribute_methods('rankdir')
    graph.create_attribute_methods('splines')
    graph.create_attribute_methods('rank')

    graph.set_rankdir("TB")
    graph.set_splines("ortho")
    graph.set_rank("same")


    graph.set_edge_defaults(constraint='true', style='invis')

    for i in range(dim[0]):
        for j in range(dim[1]-1):
            graph.add_edge(pydot.Edge(nodes[i + dim[1] * j],nodes[i + dim[1] * (j+1)]))


    for i in range(dim[0]):
        s_graph = pydot.Subgraph(rank = 'same')
        for j in range(dim[1]):
            s_graph.add_node(pydot.Node(nodes[j + dim[1] * i]))
        graph.add_subgraph(s_graph)
 

def arc_struct(graph: nx.DiGraph):

    label_dict = {node : "pe" + str(label) for node,label in zip(graph.nodes(),graph.nodes())}
    
    # set all nodes attributes
    nx.set_node_attributes(graph,'filled','style')
    nx.set_node_attributes(graph,'square','shape')
    nx.set_node_attributes(graph,'false','fixedsize')
    nx.set_node_attributes(graph,'0.6','width')
    nx.set_node_attributes(graph,'8','fontsize')
    nx.set_node_attributes(graph,'#FFc1c1','fillcolor')
    nx.set_node_attributes(graph,label_dict,'label')
    nx.set_node_attributes(graph,' ','tooltip')

    # set all edges attributes
    nx.set_edge_attributes(graph,'false','constraint')
    nx.set_edge_attributes(graph,'penwidth(0.1)','style')
    nx.set_edge_attributes(graph,'grey89','color')

        

    return nx.nx_pydot.to_pydot(graph) , list(graph.nodes())

