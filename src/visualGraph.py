import networkx as nx 
import pydot
import pandas as pd
import matplotlib.pyplot as plt
import altair as alt
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
    plt.show()


def print_dot(dot: pydot.Dot) -> None:



    text_file = open("Output.dot", "w")
    text_file.write(string_dot)
    text_file.close()


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
    with open("output.dot", 'w') as f:
        f.write(graph_string)
    f.close()



 

def arch_struct(graph: nx.DiGraph):

    label_dict = {node : "pe" + str(label) for node,label in zip(graph.nodes(),graph.nodes())}
    
    # set all nodes attributes
    nx.set_node_attributes(graph,'filled','style')
    nx.set_node_attributes(graph,'square','shape')
    nx.set_node_attributes(graph,'false','fixedsize')
    nx.set_node_attributes(graph,'0.6','width')
    nx.set_node_attributes(graph,'8','fontsize')
    nx.set_node_attributes(graph,'#FFc1c1','fillcolor')
    nx.set_node_attributes(graph,label_dict,'label')
    nx.set_node_attributes(graph,' ','tooltip')

    # set all edges attributes
    nx.set_edge_attributes(graph,'false','constraint')
    nx.set_edge_attributes(graph,'penwidth(0.1)','style')
    nx.set_edge_attributes(graph,'grey89','color')

        

    return nx.nx_pydot.to_pydot(graph) , list(graph.nodes())

