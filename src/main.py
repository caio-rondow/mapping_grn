
from matplotlib.pyplot import plot
from mappingGRN import mappingGRN
import networkx as nx
from grn2dot.grn2dot import Grn2dot
import json2graph as json2graph
import pandas as pd
import matplotlib.pyplot as plt
import visualGraph as visualGraph



def main():
    grn2dot = Grn2dot('misc/expressions.ALL.txt')
    GRN = grn2dot.get_nx_digraph()

    
    
    mapping = mappingGRN('misc/mesh_8x8.json', GRN)
    mapping.simulated_annealing()

    visualGraph.sa_curve(mapping.get_allcost())

    
    # dot,nodes = visualGraph.arc_struct(mapping.get_nx_arc())
    # visualGraph.build_dot(dot,nodes,[15,15])
    # visualGraph.print_dot(dot)




    


    


if __name__ == '__main__':
    main()