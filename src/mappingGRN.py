from asyncio import new_event_loop
from dbm import dumb
from opcode import opname
from tokenize import Double
from turtle import position

from numpy import double, mat
import json2graph
import networkx as nx
import json
import math 
import random as rand
from tqdm import tqdm

class mappingGRN:

    _instance = None

    def __init__(self, file_path, graph) -> None:
        self.set_cgra(file_path)
        self.set_grn(graph)



    # SETS
    def set_cgra(self, file_path) -> None:
        f               = open(file_path)
        self.cgra       = json2graph.make_digraph(json.load(f))
        self.arc_size   = self.cgra.__len__()
        nx.set_edge_attributes(self.cgra,1,'weight')
        nx.set_edge_attributes(self.cgra,'penwidth(0.1)','style')
        nx.set_edge_attributes(self.cgra,'grey89','color')
        nx.set_edge_attributes(self.cgra,' ','tooltip')


        nx.set_node_attributes(self.cgra,'8','fontsize')
        nx.set_node_attributes(self.cgra,'#FFc1c1','fillcolor')
        nx.set_node_attributes(self.cgra,' ','label')
        nx.set_node_attributes(self.cgra,' ','tooltip')

        f.close()


    def set_grn(self, graph: nx.DiGraph ) -> None:
        # init values
        self.wcase = self.cost = 0
        self.ctSwap = 0
        self.allCost=[]
        self.r_mapping = {}
        self.grn = graph


        self.__random_mapping()

    # GETS
    def get_arc_size(self) -> int:
        return self.arc_size

    def get_cgra(self) -> nx.DiGraph:
        return self.cgra


    def get_grn(self) -> nx.DiGraph:
        return self.grn


    def get_mapped_grn(self) -> dict:
        ''' Return r_mapping
        '''
        return self.r_mapping


    def get_worstcase(self) -> int:
        return self.wcase


    def get_allcost(self) -> int:
        return self.allCost

    
    def get_num_swaps(self) -> int:
        return self.ctSwap


    def get_dot(self):
        return json2graph.nx_2_dot(self.cgra)

    def display_arc(self):
        bline = math.sqrt(self.arc_size)
        for i in range(self.arc_size):
            if i%bline==0 : print()
            node = self.__arc_2_grn(i)
            if(self.grn.has_node(node)):
                print(node[1], end=' ')
            else: 
                print('-', end=' ')
        print()


    # METHODS
    def is_monomorphism(self, GRN: nx.DiGraph) -> bool:
        m = nx.algorithms.isomorphism.DiGraphMatcher(self.cgra,GRN)
        return m.subgraph_is_monomorphic()


    def get_monomorphism(self, GRN: nx.DiGraph) -> list:
        if(self.is_monomorphism(GRN)):
            m = nx.algorithms.isomorphism.DiGraphMatcher(self.cgra,GRN)
            m_list = list(m.subgraph_monomorphisms_iter())
        return m_list




    def __random_mapping(self) -> None: 
        """ Return a dictionary where the keys are nodes in the architecture and the values are random nodes from the graph.

            Parameters
            ----------
            graph: digraph
                A networkX digraph
                A gene regulatory network that will be used as values for de dictionary.

            Notes
            ----------  
        """
        # create a list with all pes
        empty_pe = []
        for i in range(self.arc_size):
            empty_pe.append(i)
        
        # choose random values [0, arcSize_nXn) to map the grn nodes in it
        arc_nodes = rand.sample(empty_pe, len( self.grn.nodes() ) )

        # map PE : NODE
        for node, k in zip(self.grn.nodes(), arc_nodes):
            self.r_mapping[k] = node


    def __grn_2_arc(self,node):
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


    def __arc_2_grn(self,int):
        try:    return self.r_mapping[int]
        except: return int

    def graph_visu(self) -> nx.DiGraph:
        """ Graph visualization with .dot 

            1. Rodar o simulated annealing
            2. pegar todas as conexões da grn e rodar elas no crga
            2.1. para cada gene -> key referente no dict (source e target)
            2.2. rodar dijkstra_path


            Notes
            ---------
            Ainda a fazer;
            1. rodar dijkstra_path e salvar cada custo não repetido em um dict sendo:
                key = custo e valule = cor

                color[distance] -> #FF0000

        """
        self.simulated_annealing()

        for node in self.grn.nodes():
            arch_node = self.__grn_2_arc(node)
            nx.set_node_attributes(self.cgra, {arch_node: {'label':node}})

                
        for node in self.grn.nodes():
            pe_source = self.__grn_2_arc(node)
            for neighbor in self.grn.neighbors(node):
                pe_target = self.__grn_2_arc(neighbor)
                distance,path = nx.single_source_dijkstra(self.cgra,pe_source,pe_target)
                color = ["#"+''.join([rand.choice('ABCDEF0123456789') for i in range(6)])]
                if distance <= self.get_worstcase() * 0.5: continue
                for path_node,i in zip(path,range(len(path))):
                    if path_node == path[-1]: break
                    self.cgra[path_node][path[i+1]]['color'] = color[0]
                    self.cgra[path_node][path[i+1]]['tooltip'] = "{} to {}".format(node,neighbor)

    
        return self.get_cgra()

        

       

    def total_edge_cost(self) -> int:
        """ Returns the init_cost edge cost from peX to peY.
            Also calculates the worst case cost. """

        # Reset costs
        self.cost,self.wcase=0,0
        for edge in self.grn.edges():
            # Get edge xy from grn
            x = self.__grn_2_arc(edge[0])
            y = self.__grn_2_arc(edge[1])

            # Calcualte distance between peX and peY
            dist_xy = nx.dijkstra_path_length(self.cgra,x,y)
            self.cost += dist_xy

            # Calculate worst case
            if dist_xy > self.wcase: self.wcase = dist_xy

        return self.cost


    def __evaluate_move(self,u,v,peU,peV) -> int:
        """ Returns the local cost from peU to all neighbors peW and
            the new local cost from peU (where peU is on peV) to
            to all neighbors peW. """

        localC,newLocalC=0,0 
        if (self.grn.has_node(u)==True):
            for w in self.grn.neighbors(u):
                if w==u: continue # Calculate distance only for the neighbors of v
                peW = self.__grn_2_arc(w)
                localC      += nx.dijkstra_path_length(self.cgra,peU,peW)
                newLocalC   += nx.dijkstra_path_length(self.cgra,peV,peW)    
            
            for w in self.grn.predecessors(u):
                if w==u: continue # Calculate distance only for the neighbors of v
                peW = self.__grn_2_arc(w)
                localC      += nx.dijkstra_path_length(self.cgra,peW,peU)
                newLocalC   += nx.dijkstra_path_length(self.cgra,peW,peV)

        return localC, newLocalC


    def __switching_cost(self,u,v,peU,peV,init_cost) -> int:
        """ Return the new cost from peU to peV """

        uLocal_cost     = 0 # -> cost from peU to all neighbors of node u
        uNew_local_cost = 0 # -> cost from peV to all neighbors of node u
        
        vLocal_cost     = 0 # -> cost from peV to all neighbors of node v
        vNew_local_cost = 0 # -> cost from peU to all neighbors of node v
        
        local_cost      = 0 # -> total local cost     (uLocal_cost + vLocal_cost)
        new_local_cost  = 0 # -> total new local cost (uNew_local_cost + vNew_local_cost)  

        # get partial local costs and new local costs
        uLocal_cost, uNew_local_cost = self.__evaluate_move(u,v,peU,peV)
        vLocal_cost, vNew_local_cost = self.__evaluate_move(v,u,peV,peU)

        # get total local cost and new local_cost
        local_cost      = uLocal_cost + vLocal_cost
        new_local_cost  = uNew_local_cost + vNew_local_cost 

        # new cost
        return (init_cost - local_cost + new_local_cost)


    def __randpes(self, inf, sup):
        # Choose a random pe between inf and sup
        peU = rand.randint(inf,sup)
        u   = self.__arc_2_grn(peU)

        # Case peU is empty
        if( self.grn.has_node(u)==False ):
            # Choose a random grn node
            v = rand.choice( list( self.r_mapping.values()) )
            # and find it in arc
            peV = self.__grn_2_arc(v)
        else:
            # Choose a random pe between inf and sup
            peV = rand.randint(inf,sup)

        return peU, peV

    def __range(self,max,dec,min) -> int:
        a_n = math.log10(min)
        a_1 = math.log10(max)
        q = math.log10(dec)

        n =  ((a_n -  a_1) + q) / q

        return int(n)
    


    def simulated_annealing(self) -> None:
        """ 
            Aplies Simulated Annealing algorithm on a GRN mapped into CGRA
            - Starts with a random mapped GRN
            - Expected to end up with a lower cost mapped GRN  

            Template
            --------
            let 'T' be the temperature, 'init_cost' the total edge cost and 'n' the arc length.

            T           <- 100
            init_cost   <- total_edge_cost()

            while T>0.00001:
                choose random pe's:
                    peU,peV <- rand( [0, n²) ), rand( [0, n²) )

                map pe's to grn nodes:
                    u,v <- CGRA_2_GRN(peU), CGRA_2_GRN(peV)

                Calculate new cost:
                    if u is a node from grn then:
                        evaluate_move(u,v,peU,peV)                
                    if v is a node from grn then:
                        evaluate_move(v,u,peV,peU)
                    new_cost <- init_cost - local_cost + new_local_cost

                Calculate acceptance probability:
                    accProb <- exp(-1 * (dC/T) )

                if new_cost < init_cost or rand([0,...,1]) < accProb then:
                    make a swap between peU and peV

                decrease temperature
        """

        # INIT
        T=100                               # Start Simulated Annealing temperature
        init_cost=self.total_edge_cost()    # Calculate current init_cost edge cost
        # interval of pe's
        inf = 0
        sup = self.arc_size-1
        n_range = self.__range(T,0.999,0.00001)

        for interation in tqdm(
            range(n_range),
            position=0,
            leave = True,
            desc= f"Simulated Annealing: "
        ):
            # Choose random Pe's
            peU, peV = self.__randpes(inf,sup)
            
            # map pe's to grn nodes
            u = self.__arc_2_grn(peU)
            v = self.__arc_2_grn(peV)

            # Calculate new cost 
            new_cost = self.__switching_cost(u,v,peU,peV,init_cost)

            # Calculate acceptance probability
            dC      = abs(new_cost - init_cost) # variation of cost, delta C
            accProb = math.exp(-1 * (dC/T) )

            # If it's a good swap
            is_good = new_cost<init_cost or rand.random()<accProb
            if is_good:
                # Swap peU content with peV content
                self.r_mapping.update({peU:v, peV:u})
                # progression of costs and num. of swaps
                if self.ctSwap%8==0: 
                    self.allCost.append([self.total_edge_cost(),self.ctSwap])
                self.ctSwap += 1

            # Decrease temp 
            T *= 0.999








# =========================== fix needed =========================== #



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