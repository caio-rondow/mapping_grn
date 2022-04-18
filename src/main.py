from mappingGRN import mappingGRN
import networkx as nx
from grn2dot.grn2dot import Grn2dot
import json2graph
from pyvis import  network as net

def main():
    grn2dot = Grn2dot('misc/Benchmark_188.txt')
    GRN = grn2dot.get_nx_digraph()
    
    mapping = mappingGRN('misc/mesh_15x15.json', GRN)
    
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


if __name__ == '__main__':
    main()