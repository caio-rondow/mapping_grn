from platform import node
from matplotlib import style
import networkx as nx
import pydot


def build_dot(graph: pydot.Dot) -> pydot.Dot:

    graph.create_attribute_methods('rankdir')
    graph.create_attribute_methods('splines')

    graph.set_rankdir("TB")
    graph.set_splines("ortho")
    graph.set_node_defaults(style='filled', shape='square', fixedsize='false', width='0.6')

    for node,i in zip(graph.get_nodes(),range(len(graph.get_nodes()))):
        graph.get_node(node)[0].set_label(i)
        graph.get_node(node)[0].set_tooltip(" ")
        graph.get_node(node)[0].set_fontsize("8")
        graph.get_node(node)[0].set_fillcolor("#FFc1c1")

    graph.set_edge_defaults(constraint="false", style="penwidth(0.1)", color="grey89")


def modelando_grafo(graph: nx.DiGraph) -> None:
    

    nx.set_node_attributes(graph,'filled','style')
    nx.set_node_attributes(graph)


