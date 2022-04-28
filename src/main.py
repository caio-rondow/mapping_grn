from mappingGRN import mappingGRN
from grn2dot.grn2dot import Grn2dot
import visualGraph as visualGraph




def main():
    grn2dot = Grn2dot('misc/expressions.ALL.txt')
    GRN = grn2dot.get_nx_digraph()

    mapping = mappingGRN('misc/mesh_8x8.json', GRN)
    #mapping.simulated_annealing()

    ### GRAPH TOTAL_COST X N_SWAPS ###
    # visualGraph.sa_curve(mapping.get_allcost())
    
    ### CGRA VISUALIZATION ###
    dot,nodes = visualGraph.arch_struct(mapping.get_cgra())
    visualGraph.build_dot(dot,nodes,[8,8])


    ### HISTOGRAM OF N TIMES A PE WAS USED ###
    # visualGraph.num_pes_used(10,mapping,GRN)


if __name__ == '__main__':
    main()