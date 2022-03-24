from mappingGRN import mappingGRN
import networkx as nx
from grn2dot.grn2dot import Grn2dot

def adjust_node(G: nx.DiGraph,node,out_degree,i) -> None:
    """ Create auxiliary nodes in a GRN if the out_degree of a gene is up to 4, 
        connecting this gene with the new node

        Parameters
        ----------
        GRN: digraph
            A networkX digraph
            A gene regulatory network that will be modeled into a new graph with max out_degree = 4

        Returns
        ----------
        The GRN with new nodes used as bridges for high out_degree genes, each new node have a atribute
        'bridge' = True

        Notes
        ----------  
    """
    
    if(G.degree(node) <= 4 or out_degree == 0):
        return
    
    G.add_node(i,bridge = True)
    G.add_edge(node, i)

    remove_edges = []
    degree = G.degree(node)
    for neighbor in G.neighbors(node):
        if degree <= 4: break
        G.add_edge(i,neighbor)
        remove_edges.append([node,neighbor])
        degree -= 1
        
    G.remove_edges_from(remove_edges)

    adjust_node(G,i,G.out_degree(i),i+1)





def adjust_GRN(GRN: nx.DiGraph) -> nx.DiGraph:
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
    nodes = GRN.nodes()

    for node,i in zip(nodes,range(len(nodes))):
        if(G.degree(node) > 4): adjust_node(G,node,G.out_degree(node),i+11)
    
    return G

def funcao_custo(graph):

def main():
    grn2dot = Grn2dot('misc\expressionsALL.txt')
    mapping = mappingGRN('misc\ex.json')
    rato = grn2dot.get_nx_digraph()


    grid = nx.grid_2d_graph(8,8)
    grid = nx.to_directed(grid)
    rato = grn2dot.digraph()

    matcher = nx.algorithms.isomorphism.DiGraphMatcher(grid,rato)
    b = matcher.subgraph_is_monomorphic()
    print(b)
    if(b): 
        c = matcher.subgraph_monomorphisms_iter()

    print(list(c))

if __name__ == '__main__':
    main()