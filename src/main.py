from matplotlib.pyplot import plot
from numpy import dtype
from mappingGRN import mappingGRN
import networkx as nx
from grn2dot.grn2dot import Grn2dot
import json2graph as json2graph
from pyvis import  network as net
import pandas as pd
import matplotlib.pyplot as plt
import pydot




def main():
    grn2dot = Grn2dot('misc/Benchmark_5.txt')
    GRN = grn2dot.get_nx_digraph()
    
    mapping = mappingGRN('misc/mesh_8x8.json', GRN)
    #mapping.simulated_annealing()

    # data = mapping.get_allcost()

    # df = pd.DataFrame(data)

    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # ax.plot(df[1], df[0], color='tab:blue')

    # plt.show()

    arc = mapping.get_nx_arc()



    


if __name__ == '__main__':
    main()