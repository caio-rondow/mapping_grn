from mappingGRN import mappingGRN
import networkx as nx
from grn2dot.grn2dot import Grn2dot
import json2graph
from pyvis import  network as net

def main():
    grn2dot = Grn2dot('misc/Benchmark_5.txt')
    GRN = grn2dot.get_nx_digraph()
    mapping = mappingGRN('misc/mesh_8x8.json', GRN)
    #G = mapping.adjust_GRN(GRN)
    dict = mapping.random_mapping()
    print(dict)
    mapping.simulated_annealing()
    wcase   = mapping.get_worstcase()
    cost    = mapping.total_edge_cost()
    print(dict)
    print('cost ::{}'.format(cost))
    print('wcase::{}'.format(wcase))


if __name__ == '__main__':
    main()