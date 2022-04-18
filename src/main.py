from mappingGRN import mappingGRN
from visualGraph import sa_curve
from grn2dot.grn2dot import Grn2dot
from pyvis import  network as net
import networkx as nx
import json2graph

def printMap(mapping):
    print('==================================')
    print('INITIAL ARC', end='')
    mapping.display_arc()

    Co = mapping.total_edge_cost()
    Wo = mapping.get_worstcase()

    print( 'Co::{}'.format(Co) )
    print( 'Wo::{}'.format(Wo) )
    
    # Applying SA
    mapping.simulated_annealing()

    print('\nAPPYLING SA...', end='')
    mapping.display_arc()

    Cf = mapping.total_edge_cost()
    Wf = mapping.get_worstcase()

    print( 'Cf::{}'.format(Cf) )
    print( 'Wf::{}'.format(Wf) )
    
    dC = abs(Cf-Co)
    dW = abs(Wf-Wo)
    nS = mapping.get_num_swaps()
    aC = mapping.get_allcost()

    print( '\ndC::{}'.format(dC) )
    print( 'dW::{}'.format(dW) )

    print( 'nS::{}'.format(nS) )
    print( 'aC::{}'.format(aC) )
    
    print('==================================') 

def main():
    grn2dot = Grn2dot('misc/Benchmark_5.txt') 
    GRN = grn2dot.get_nx_digraph()
    
    mapping = mappingGRN('misc/mesh_8x8.json', GRN)
    #printMap(mapping)

    mapping.simulated_annealing()
    data = mapping.get_allcost()
    sa_curve(data)

if __name__ == '__main__':
    main()