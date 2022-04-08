from asyncio import new_event_loop
from opcode import opname
from turtle import position
import json2graph
import networkx as nx
import json
from pyvis import network as net
import math 
import random as rand

class mappingGRN:

    _instance = None

    def __init__(self,file_path,graph) -> None:
        f               = open(file_path)
        self.g1         = json2graph.make_digraph(json.load(f))
        self.r_mapping  = {}
        self.grn        = graph
        self.wcase      = 0
        self.cost       = 0
        
        nx.set_edge_attributes(self.g1,1,'weight')

    def get_digraph(self) -> nx.DiGraph:
        return self.digraph
    
    def set_grn(self,graph: nx.DiGraph()) -> None:
        self.grn = graph
        self.r_mapping = {}
        self.random_mapping()
        self.wcase,self.cost=0,0

    def get_dot(self) -> None:
        print(json2graph.nx_2_dot(self.g1))

    # Graph visualization with Pyvis
    def graph_visu(self,__directed=True) -> None:
        # Make a .html for download (nt.show() has some kind of problem on Colab)
        nt = net.Network(directed = __directed)
        nt.from_nx(self.digraph)
        nt.save_graph('digraph_visualization.html')
        # Download the .html
        nt.show('digraph_visualization.html')

    def adjust_node(self,G: nx.DiGraph,node,out_degree,i=0) -> None:
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
        m = nx.algorithms.isomorphism.DiGraphMatcher(self.g1,GRN)
        return m.subgraph_is_monomorphic()

    def get_monomorphism(self,GRN: nx.DiGraph) -> list:
        if(self.is_monomorphism(GRN)):
            m = nx.algorithms.isomorphism.DiGraphMatcher(self.g1,GRN)
            m_list = list(m.subgraph_monomorphisms_iter())
        
        return m_list

    def random_mapping(self) -> dict:
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
        for arc_node,grn_node in zip(self.g1.nodes(),self.grn.nodes()):
            self.r_mapping[arc_node] = grn_node
        return self.r_mapping 

    def grn_2_arc(self,node):
        """ Give one node in the GRN, return the CGRA's node that it is in.

            Parameters
            ----------
            node: string
                A node in the GRN graph

            Returns
            ----------
            A node in the CGRA.

            Notes
            -
        """
        key_list = list(self.r_mapping.keys())
        val_list = list(self.r_mapping.values())

        try:    position = val_list.index(node)
        except: return node
        return key_list[position]

    def arc_2_grn(self,int):
        try:    return self.r_mapping[int]
        except: return int
       
    def get_worstcase(self):
        return self.wcase

    def total_edge_cost(self) -> int:
        """
            Returns the total edge cost from peX to peY.
            Also calculates the worst case cost.
        """
        # Reset costs
        self.cost,self.wcase=0,0
        for edge in self.grn.edges():
            # Get edge xy from grn
            x = self.grn_2_arc(edge[0])
            y = self.grn_2_arc(edge[1])

            # Calcualte distance between peX and peY
            dist_xy = nx.dijkstra_path_length(self.g1,x,y)
            self.cost += dist_xy

            # Calculate worst case
            if dist_xy > self.wcase: self.wcase = dist_xy
        return self.cost

    def evaluate_move(self,u,v,peU,peV) -> int:
        """
            Returns the local cost from peU to all neighbors peW and
            the new local cost from peU (where peU is on peV) to
            to all neighbors peW.
        """
        localC,newLocalC=0,0 
        if (self.grn.has_node(u)==True):
            for w in self.grn.neighbors(u):
                if w==u: continue # Calculate distance only btw neighbors of v
                peW = self.grn_2_arc(w)
                localC      += nx.dijkstra_path_length(self.g1,peU,peW)
                newLocalC   += nx.dijkstra_path_length(self.g1,peV,peW)    
            for w in self.grn.predecessors(u):
                if w==u: continue # Calculate distance only btw neighbors of v
                peW = self.grn_2_arc(w)
                localC      += nx.dijkstra_path_length(self.g1,peW,peU)
                newLocalC   += nx.dijkstra_path_length(self.g1,peW,peV)
        return localC, newLocalC

    def simulated_annealing(self) -> None:
        """ 
            Aplies Simulated Annealing algorithm on a GRN mapped into CGRA
            - Starts with a random mapped GRN
            - Expected to end up with a lower cost mapped GRN  
        """
        T=100                        # Start Simulated Annealing temperature
        total=self.total_edge_cost() # Calculate current total edge cost

        while(T>0.00001):
            # Choose a random Pe between [0, mesh_nXn-1]
            peU,peV = rand.randint(0, self.g1.__len__()-1 ), rand.randint(0, self.g1.__len__()-1 )

            # Get a node from GRN mapped on CGRN
            u,v = self.arc_2_grn(peU), self.arc_2_grn(peV)

            # Calculate new Cost newC
            lc1,lc2,nlc1,nlc2=0,0,0,0
            if (self.grn.has_node(u)==True):
                lc1,nlc1 = self.evaluate_move(u,v,peU,peV)
            if (self.grn.has_node(v)==True):
                lc2,nlc2 = self.evaluate_move(v,u,peV,peU)
            lc      = lc1+lc2
            nlc     = nlc1+nlc2 
            newC    = total-lc+nlc      
            
            # Calculate acceptance probability
            dE      = abs(newC - total)
            accProb = math.exp(-1 * (dE/T) )
            
            # If it's a good swap
            if(newC < total or rand.random() < accProb):
                # Swap peU content with peV content
                self.r_mapping.update({peU:v, peV:u})

            # Decrease temp 
            T *= 0.999