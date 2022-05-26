from grn2dot.grn2dot import Grn2dot
import src.include.visualGraph as visualGraph
from src.mappingGRN import mappingGRN
import src.algorithms.simulated_anealling as sm
from src.include.save_script import save_script
import src.algorithms.simulated_anealling2t as sm2t





def main():
    grn2dot = Grn2dot('misc/Benchmark_53.txt')

    GRN = grn2dot.get_nx_digraph()

    arch_path = 'misc/arch/15x15/cgra_mesh_ho_15x15.json'

    aux = arch_path.split('/')
    aux = aux[1].split('.')
    fname = aux[0]


    mapping = mappingGRN(arch_path, GRN)
    # mapping2t = mappingGRN(arch_path, GRN)
    
    # print(mapping.get_mapped_grn())
    # print(mapping2t.get_mapped_grn())

    # print("SM")
    # sm.simulated_annealing(mapping,data=True)
    # print("SM2T")
    # sm2t.simulated_annealing(mapping,data=True)


    # ### CGRA VISUALIZATION ###
    arch = mapping.graph_visu()
    dot,nodes = visualGraph.arch_struct(arch)
    visualGraph.build_dot(dot,nodes,[15,15],fname)

    # ### GRAPH TOTAL_COST X N_SWAPS ###
    # visualGraph.sa_curve(mapping.get_allcost(),fname)

    # # ### HISTOGRAM OF N TIMES A PE WAS USED ###
    # # # visualGraph.num_pes_used(10,mapping,GRN)

    # # print(mapping.get_num_swaps())

    # save_script("misc\\grn_benchmarks-main","misc\\arch\\15x15")





if __name__ == '__main__':
    main()