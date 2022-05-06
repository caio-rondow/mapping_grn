from mappingGRN import mappingGRN
from grn2dot.grn2dot import Grn2dot
import visualGraph as visualGraph




def main():
    grn2dot = Grn2dot('misc/Benchmark_53.txt')
    GRN = grn2dot.get_nx_digraph()

    mapping = mappingGRN('misc/mesh_8x8_het.json', GRN)
    # mapping.simulated_annealing()

    

    ### CGRA VISUALIZATION ###
    arch = mapping.graph_visu()
    dot,nodes = visualGraph.arch_struct(arch)
    visualGraph.build_dot(dot,nodes,[8,8])

    ### GRAPH TOTAL_COST X N_SWAPS ###
    visualGraph.sa_curve(mapping.get_allcost())

    # ### HISTOGRAM OF N TIMES A PE WAS USED ###
    # # visualGraph.num_pes_used(10,mapping,GRN)

    # print(mapping.get_num_swaps())




if __name__ == '__main__':
    main()