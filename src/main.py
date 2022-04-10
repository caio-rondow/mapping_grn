from mappingGRN import mappingGRN
import networkx as nx
from grn2dot.grn2dot import Grn2dot
import json2graph
from pyvis import  network as net

def main():
    grn2dot = Grn2dot('misc/B_bronchiseptica_and_T_retortaeformis_coinfectio_53.txt')

    GRN = grn2dot.get_nx_digraph()
    
    mapping = mappingGRN('misc/mesh_8x8.json', GRN)
    
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

    print( '\ndC::{}'.format(dC) )
    print( 'dW::{}'.format(dW) )
    print('==================================')

if __name__ == '__main__':
    main()