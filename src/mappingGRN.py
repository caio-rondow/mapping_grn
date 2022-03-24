from opcode import opname
import json2graph
import networkx as nx
import json
from pyvis import  network as net


class mappingGRN:

    _instance = None

    def __init__(self,file_path) -> None:
        f = open(file_path)
        self.digraph = json2graph.make_digraph(json.laod(f))


    # Graph visualization with Pyvis
    def graph_visu(self,__directed=True) -> None:
        # Make a .html for download (nt.show() has some kind of problem on Colab)
        nt = net.Network(directed = __directed)
        nt.from_nx(self.digraph)
        nt.save_graph('digraph_visualization.html')
        # Download the .html
        nt.show('digraph_visualization.html')


    # def adjust_GRN(GRN: nx.DiGraph()) -> None:
    #     """ Create auxiliary nodes in a GRN if the out_degree of a gene is up to 4, 
    #         connecting this gene with the new node

    #         Parameters
    #         ----------
    #         GRN: digraph
    #             A networkX digraph
    #             A gene regulatory network that will be modeled into a new graph with max out_degree = 4

    #         Returns
    #         ----------
    #         The GRN with new nodes used as bridges for high out_degree genes, each new node have a atribute
    #         'bridge' = True

    #         Notes
    #         ----------  
    #     """
        
    #     nodes = GRN.nodes()

    #     for node in nodes:
    #         node_out_degree = GRN.out_degree(node)
    #         if(node_out_degree < 4):
    #             new_node_out_degree = node_out_degree - 4
    #             GRN.add_node(-1*node, bridge = True)
    #             GRN.add_edges(-1*node,node)

            










