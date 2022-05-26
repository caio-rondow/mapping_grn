

from numpy import empty
import src.include.json2graph as json2graph
import src.algorithms.simulated_anealling as sa
import networkx as nx
import json
import math 
import random as rand

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
        nx.set_node_attributes(self.cgra,'#FFFFFF','fillcolor')
        nx.set_node_attributes(self.cgra,' ','label')
        nx.set_node_attributes(self.cgra,' ','tooltip')
        nx.set_node_attributes(self.cgra,'square','shape')

        f.close()


    def set_grn(self, graph: nx.DiGraph ) -> None:
        # init values
        self.wcase = self.cost = 0
        self.ctSwap = 0
        self.allCost=[]
        self.r_mapping = {}
        self.grn = graph


        self.__random_mapping(4)

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
            node = self.arc_2_grn(i)
            if(self.grn.has_node(node)):
                print(node[1], end=' ')
            else: 
                print('-', end=' ')
        print()

    def get_all_stats(self) -> None:
        print(
            f"{'Number of PEs in CGRA:' : <30}{self.get_arc_size() : >10}",
            f"\n{'Number of genes in the GRN:' : <30}{self.get_grn().number_of_nodes() : >10}",
            f"\n{'Total number of swaps:' : <30}{self.get_num_swaps() : >10}",
            f"\n{'Total cost:' : <30}{self.total_edge_cost() : >10}",
            f"\n{'Worst path cost:' : <30}{self.get_worstcase() : >10}"
        )
                                                                                

    # METHODS
    def is_monomorphism(self, GRN: nx.DiGraph) -> bool:
        m = nx.algorithms.isomorphism.DiGraphMatcher(self.cgra,GRN)
        return m.subgraph_is_monomorphic()


    def get_monomorphism(self, GRN: nx.DiGraph) -> list:
        if(self.is_monomorphism(GRN)):
            m = nx.algorithms.isomorphism.DiGraphMatcher(self.cgra,GRN)
            m_list = list(m.subgraph_monomorphisms_iter())
        return m_list


    def __random_mapping(self, seed=None) -> None: 
        """ Return a dictionary where the keys are nodes in the architecture and the values are random nodes from the graph.

            Parameters
            ----------
            graph: digraph
                A networkX digraph
                A gene regulatory network that will be used as values for de dictionary.

            Notes
            ----------  
        """
        if seed != None:
            rand.seed(seed)

        empty_pe = []
        # create a list with all pes

        if self.grn.number_of_nodes() < 64:
            for i in range(2,10):
                for j in range(2,10):
                    empty_pe.append(15 * i + j) # fixado para arch 15x15
            print(empty_pe)
        else:
            empty_pe = list(range(self.arc_size))
        
        # choose random values [0, arcSize_nXn) to map the grn nodes in it
        arc_nodes = rand.sample(empty_pe, len( self.grn.nodes() ) )

        # map PE : NODE
        for node, k in zip(self.grn.nodes(), arc_nodes):
            self.r_mapping[k] = node


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


    def arc_2_grn(self,node):
        try: 
            return self.r_mapping[node]   
        except: return None
 
        

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
        sa.simulated_annealing(self,data=True)

        # list of lists where each one is a path in the GRN on CGRA
        paths = []
        colors = {}
        

        # dict that count 0. how many times an PE is source
        #                 1. how many times an PE is target
        #                 2. how many times an PE is pasage for each PE
        # dict[PE] = [0,1,2]
        pe_stats = {pe : [0,0,0] for pe in range(self.get_arc_size())}
                
        for node in self.grn.nodes():
            pe_source = self.grn_2_arc(node)
            pe_stats[pe_source][0] += 1
            for neighbor in self.grn.neighbors(node):
                pe_target = self.grn_2_arc(neighbor)
                pe_stats[pe_target][1] += 1
                distance,path = nx.single_source_dijkstra(self.cgra,pe_source,pe_target)
                if distance <= self.get_worstcase() * 0.5: continue
                if distance not in colors.keys():
                    while True:
                        color = ["#"+''.join([rand.choice('ABCDEF0123456789') for i in range(6)])]
                        if color not in colors.values(): break
                    colors[distance] = color
                paths.append([path,distance])
            nx.set_node_attributes(self.cgra,{pe_source: {'fillcolor':'#666666'}}) #fill color for PEs with a gene in
            

        # add edge attributes:
        #                     a) colors by path's distance
        #                     b) tooltip showing path's source and target: pe_name(gene_name) to pe_name(gene_name) 
        for path,distance in paths:
            for path_node,i in zip(path,range(len(path))):
                if path_node == path[-1]: break
                if i > 0: pe_stats[path_node][2] += 1
                self.cgra[path_node][path[i+1]]['color'] = colors[distance][0]
                self.cgra[path_node][path[i+1]]['tooltip'] = "{}({}) to {}({})".format(path[0],
                                                                                        self.arc_2_grn(path[0]),
                                                                                        path[-1],
                                                                                        self.arc_2_grn(path[-1]))

        # add node attributes:
        #                     a) fillcolor of PEs used as bridge
        #                     b) label showing PEs stats: in degree, out degree, bridge count
        #                     c) tooltip showing gene's name in PE
        for pe in self.cgra.nodes():
            grn_node = self.arc_2_grn(pe)
            if (pe_stats[pe][0] == 0 and pe_stats[pe][1] == 0 and pe_stats[pe][2] != 0):
                nx.set_node_attributes(self.cgra,{pe: {'fillcolor':'#bdbdbd'}}) #fill color for PEs used as bridge
            nx.set_node_attributes(self.cgra, {pe: {'label':pe_stats[pe]}})
            nx.set_node_attributes(self.cgra,{pe: {'tooltip':f"name: {grn_node},\nin_degree: {self.cgra.in_degree(pe)},\nout_degree: {self.cgra.out_degree(pe)}" }})


        bridge_dict = {}
        for k, v in pe_stats.items():
            bridge_dict[k] = (v[0] + v[1] + v[2])

        bridge_dict = {k: v for k, v in sorted(bridge_dict.items(), key=lambda item: item[1])}
        values = list(bridge_dict.values())
        max = values[-1]

        for pe in bridge_dict.keys():
            if bridge_dict[pe] == max : 
                nx.set_node_attributes(self.cgra,{pe: {'shape': 'Msquare'}})



    
        return self.get_cgra()

        

       

    def total_edge_cost(self) -> int:
        """ Returns the init_cost edge cost from peX to peY.
            Also calculates the worst case cost. """

        # Reset costs
        self.cost,self.wcase=0,0
        for edge in self.grn.edges():
            # Get edge xy from grn
            x = self.grn_2_arc(edge[0])
            y = self.grn_2_arc(edge[1])

            # Calcualte distance between peX and peY
            dist_xy = nx.dijkstra_path_length(self.cgra,x,y)
            if dist_xy == 3 or dist_xy == 2:
                self.cost += 1
            else:
                self.cost += dist_xy

            # Calculate worst case
            if dist_xy > self.wcase: self.wcase = dist_xy

        return self.cost

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