from mappingGRN import mappingGRN
import networkx as nx
from grn2dot.grn2dot import Grn2dot
import json2graph
from pyvis import  network as net

def main():
    grn2dot = Grn2dot('misc/Benchmark_5.txt')
    mapping = mappingGRN('misc/mesh_8x8.json')
    GRN = grn2dot.get_nx_digraph()
    G = mapping.adjust_GRN(GRN)
    dict = mapping.randon_mapping(GRN)

    #print(grn2dot.get_nodes_vector())
    #print(grn2dot.get_edges_dict())

    cost = mapping.total_edge_cost(GRN)
    wcase = mapping.get_worstcase()

    print(cost)
    print(wcase)
    print(dict)
    mapping.simulated_annealing(GRN)
    print(dict)


if __name__ == '__main__':
    main()