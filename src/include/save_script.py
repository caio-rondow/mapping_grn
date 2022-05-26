import src.algorithms.simulated_anealling as sa
from src.mappingGRN import mappingGRN as mapping
from grn2dot.grn2dot import Grn2dot
from glob import glob
import numpy as np
import pandas as pd
import pathlib



def GRN_paths(path) -> list:
    ''' Return a list with all nx.Digraph object of all GRN found in path
    '''

    GRN = []
    p = pathlib.Path(path)
    paths = list(p.glob('**/expressions.ALL.txt'))
    names = get_grn_names(paths)
    for path in paths:
        grn2dot = Grn2dot(path)
        GRN.append(grn2dot.get_nx_digraph())

    return GRN,names

def get_arch_names(ARCH: list) -> list:

    arch_names = []

    for arch in ARCH:
        aux = arch.split('\\')
        aux = aux[2].split('.')
        arch_names.append(aux[0])
    
    return arch_names

def get_grn_names(GRN: list) -> list:

    grn_names = []

    for grn in GRN:
        rato = str(grn)
        aux = rato.split('\\')
        grn_names.append(aux[3])

    return grn_names

def get_data(row: list, mp: mapping) -> list:
    cost = mp.total_edge_cost()
    wc = mp.get_worstcase()

    row.extend([wc,cost])

    return row


def write_file(data: pd.DataFrame()) -> None:
    ''' Write a xls file with all data
    '''

def save_script(grn_path, arch_path):

    GRN,grn_names = GRN_paths(grn_path)
    ARCH = glob(arch_path + "/*.json")

    rows = []

    for grn in GRN:
        row = []
        if grn.number_of_nodes() > 225:
            row.extend(['-'] * (len(ARCH)*2))
            print(row)
            rows.append(row)
            continue
        for arch in ARCH:
            mp = mapping(arch,grn)
            try:
                sa.simulated_annealing(mp)
                row = get_data(row,mp)
            except:
                row.extend(['-','-'])
        print(row)
        rows.append(row)


    arch_names = get_arch_names(ARCH)
    columns = pd.MultiIndex.from_product([arch_names,['WC','Cost']])
    data = pd.DataFrame(rows,columns=columns, index=grn_names)


    try: 
        data.to_excel("output.xlsx")
    except:
        print(data)
        




