import os
import networkx as nx 
import pydot
import pandas as pd
import matplotlib.pyplot as plt
import altair as alt
from sqlalchemy import Constraint
from mappingGRN import mappingGRN



def num_pes_used(n_simu: int, mapping: mappingGRN, GRN: nx.DiGraph) -> pd.DataFrame:

    arc = mapping.get_cgra()
    data = {pe : 0 for pe in arc.nodes()}

    for i in range(n_simu):
        mapping.simulated_annealing()
        pe_dict = mapping.get_mapped_grn()
        for pe in list(pe_dict.keys()):
            if pe_dict[pe] in GRN.nodes():
                data[pe] += 1

    df = pd.DataFrame(list(data.items()), columns=['PE','N times used'])
    
    chart = alt.Chart(df).mark_bar().encode(
    alt.X("PE:N"),
    alt.Y("N times used"),
    )

    chart.save('histogram_{}PEs_{}T_{}GRN.html'.format(mapping.get_arc_size(),n_simu,GRN.number_of_nodes()))

    return df



def sa_curve(data):
    df = pd.DataFrame(data)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(df[1], df[0], color='tab:blue')
    ax.set_title('Simulated Anealling Cost')
    ax.set_xlabel('num swaps')
    ax.set_ylabel('total cost')


    path = "C:\\Users\\Windows 10\\Google Drive\\Iniciação\\codes\\mapping_grn\\benchmarks"
    file_name = f"cost_graph.svg"

    completeName = os.path.join(path, file_name)


    plt.savefig(completeName,dpi=150)
    plt.show()



def build_dot(graph: pydot.Dot, nodes: list, dim: list) -> None:

    '''
        dim = row x col -> dim[0] = row din[1] = col

        nodes[i + dim[1] * j]

        nodes[j + dim[1] * i]
    '''

    graph.create_attribute_methods('rankdir')
    graph.create_attribute_methods('splines')

    graph.set_rankdir("TB")
    graph.set_splines("ortho")

    graph_string = graph.to_string()
    graph_string = graph_string.replace("strict digraph  {", "digraph layout  {")


    # colunas
    graph_string = graph_string.replace("}","edge [constraint=true, style=invis];\n")

    for i in range(dim[0]):
        for j in range(dim[1]):
            if (j + 1) % dim[1] == 0:
                graph_string = graph_string + "%d;\n" % (nodes[i + dim[1] * j])
            else:
                graph_string = graph_string + "{} -> ".format(nodes[i + dim[1] * j])

    # linhas
    for i in range(dim[0]):
        graph_string = graph_string + "rank = same {"
        for j in range(dim[1]):
            if (j + 1) % dim[1] == 0:
                graph_string = graph_string + "%d};\n" % (nodes[j + dim[1] * i])
            else:
                graph_string = graph_string + "{} -> ".format(nodes[j + dim[1] * i])

    graph_string = graph_string + "}"

    path = "C:\\Users\\Windows 10\\Google Drive\\Iniciação\\codes\\mapping_grn\\benchmarks"
    file_name = f"{len(nodes)}nodes_{dim[0]}x{dim[1]}cgra.dot"

    completeName = os.path.join(path, file_name)

    with open(completeName, 'w') as f:
        f.write(graph_string)
    f.close()


def to_pydot(N):
    """Returns a pydot graph from a NetworkX graph N.

    Parameters
    ----------
    N : NetworkX graph
      A graph created with NetworkX

    Examples
    --------
    >>> K5 = nx.complete_graph(5)
    >>> P = nx.nx_pydot.to_pydot(K5)

    Notes
    -----
    This function is derived from NetworkX method (see more on: https://networkx.org/documentation/stable/reference/generated/networkx.drawing.nx_pydot.to_pydot.html)

    """
    import pydot

    # set Graphviz graph type
    if N.is_directed():
        graph_type = "digraph"
    else:
        graph_type = "graph"
    strict = nx.number_of_selfloops(N) == 0 and not N.is_multigraph()

    name = N.name
    graph_defaults = N.graph.get("graph", {})
    if name == "":
        P = pydot.Dot("", graph_type=graph_type, strict=strict, **graph_defaults)
    else:
        P = pydot.Dot(
            f'"{name}"', graph_type=graph_type, strict=strict, **graph_defaults
        )
    try:
        P.set_node_defaults(style = 'filled',
                            fixedsize = 'false',
                            width = '0.6')
    except KeyError:
        pass
    try:
        P.set_edge_defaults(constraint = 'false')
    except KeyError:
        pass

    for n, nodedata in N.nodes(data=True):
        str_nodedata = {k: str(v) for k, v in nodedata.items()}
        p = pydot.Node(str(n), **str_nodedata)
        P.add_node(p)

    if N.is_multigraph():
        for u, v, key, edgedata in N.edges(data=True, keys=True):
            str_edgedata = {k: str(v) for k, v in edgedata.items() if k != "key"}
            edge = pydot.Edge(str(u), str(v), key=str(key), **str_edgedata)
            P.add_edge(edge)

    else:
        for u, v, edgedata in N.edges(data=True):
            str_edgedata = {k: str(v) for k, v in edgedata.items()}
            edge = pydot.Edge(str(u), str(v), **str_edgedata)
            P.add_edge(edge)
    return P
 

def arch_struct(graph: nx.DiGraph):

    # set all nodes attributes
    # nx.set_node_attributes(graph,'square','shape')
    return to_pydot(graph) , list(graph.nodes())

