from asyncio import new_event_loop
from opcode import opname
from turtle import position
import json2graph
import networkx as nx
import json
from pyvis import  network as net


class mappingGRN:

    _instance = None

    def __init__(self,file_path) -> None:
        f = open(file_path)
        self.arc = json2graph.make_digraph(json.load(f))
        nx.set_edge_attributes(self.arc,1,'weight')


        self.r_mapping = {}
        self.cost = 0


    def get_digraph(self) -> nx.DiGraph:
        return self.digraph

    def get_dot(self) -> None:
        print(json2graph.nx_2_dot(self.arc))

    # Graph visualization with Pyvis
    def graph_visu(self,__directed=True) -> None:
        # Make a .html for download (nt.show() has some kind of problem on Colab)
        nt = net.Network(directed = __directed)
        nt.from_nx(self.digraph)
        nt.save_graph('digraph_visualization.html')
        # Download the .html
        nt.show('digraph_visualization.html')


    def adjust_node(self,G: nx.DiGraph,node,out_degree) -> None:
        """ Create auxiliary nodes in a GRN if the out_degree of a gene is up to 4, 
            connecting this gene with the new node

            Parameters
            ----------
            G: digraph
                A networkX digraph in which will be add new nodes as bridges
            node: node-like from networkx 
                a node from the graph G that will be the start link of the bride i
            out_degree: int
                the number of edges that get out of node
            i: int
                the label that will be the name of the bridge

            Returns
            ----------
            The GRN with new nodes used as bridges for high out_degree genes, each new node have a atribute
            'bridge' = True

            Notes
            ----------
            For more: https://excalidraw.com/#json=_u3d67AErZB19oQ0qeDGr,th6yvk8UgBEsTCsW0mqQKw  
        """
        
        if G.degree(node) <= 4 | out_degree == 0:
            return
        
        new_node_name = node + "_bridge"
        G.add_node(new_node_name,bridge = True)
        G.add_edge(node, new_node_name)

        remove_edges = []
        for neighbor in G.neighbors(node):
            if G.degree(node) <= 4 | out_degree == 0: break
            G.add_edge(new_node_name,neighbor)
            remove_edges.append([node,neighbor])
            
        G.remove_edges_from(remove_edges)

        self.adjust_node(G,new_node_name,G.out_degree(new_node_name))


    def merge_nodes(self,G: nx.DiGraph,merge_nodes: list) -> None:



        for m_node in merge_nodes[1:]:
            nx.contracted_nodes(G,merge_nodes[0],m_node,copy=False)
        



    def adjust_GRN(self,GRN: nx.DiGraph) -> nx.DiGraph:
        """ Create auxiliary nodes in a GRN if the out_degree of a gene is up to 4, 
            connecting this gene with the new node

            Parameters
            ----------
            GRN: digraph
                A networkX digraph
                A gene regulatory network that will be modeled into a new graph with max out_degree = 4

            Returns
            ----------
            A new networkx digraph with new nodes as bridges

            Notes
            ----------  
        """
        
        G = nx.DiGraph(GRN)

        G.remove_edges_from(nx.selfloop_edges(G))
        
        '''Fix needed'''
        
        return G
        
    def is_monomorphism(self,GRN: nx.DiGraph) -> bool:
        m = nx.algorithms.isomorphism.DiGraphMatcher(self.arc,GRN)
        return m.subgraph_is_monomorphic()

    def get_monomorphism(self,GRN: nx.DiGraph) -> list:
        if(self.is_monomorphism(GRN)):
            m = nx.algorithms.isomorphism.DiGraphMatcher(self.arc,GRN)
            m_list = list(m.subgraph_monomorphisms_iter())
        
        return m_list

    def randon_mapping(self,graph: nx.DiGraph) -> dict:
        """ Return a dictionary where the keys are nodes in the architecture and the values are random nodes from the graph.

            Parameters
            ----------
            graph: digraph
                A networkX digraph
                A gene regulatory network that will be used as values for de dictionary.

            Returns
            ----------
            A python dictionary.

            Notes
            ----------  
        """

        for arc_node,graph_node in zip(self.arc.nodes(),graph.nodes()):
            self.r_mapping[arc_node] = graph_node
    
        return self.r_mapping 


    def grn_2_arc(self,node):
        """ Give one node in the GRN, return the CGRA's node that it is in.

            Parameters
            ----------
            node: node label
                A node in the GRN graph
                Nodes can be, for example, strings or numbers. Nodes must be hashable (and not None) Python objects.

            Returns
            ----------
            key_list[position]: node label
                A node in the CGRA. If the node used as parameter is not in the GRN, returns the node itself.

            Notes
            -
        """
        key_list = list(self.r_mapping.keys())
        val_list = list(self.r_mapping.values())

        try:
            position = val_list.index(node)
        except KeyError:
            position = node

        return key_list[position]

    def arc_2_grn(self,int):
        """ Give one node in the CGRA, return the GRN's node thats is mapping in int.

            Parameters
            ----------
            int: node label
                A node in the CGRA graph
                Nodes can be, for example, strings or numbers. Nodes must be hashable (and not None) Python objects.

            Returns
            ----------
            grn_node: node label
                A node in the GRN. If the node used as parameter is not one that a GRN's node was mapped into, retunrs
                int itsef

            Notes
            -
        """
        try:
            grn_node = self.r_mapping[int]
        except KeyError:
            grn_node = int
        
        return grn_node

    def total_edge_cost(self,graph) -> int:
        for edge in graph.edges():
            x = self.grn_2_arc(edge[0])
            y = self.grn_2_arc(edge[1])

            self.cost += nx.dijkstra_path_length(self.arc,x,y)

        return self.cost

