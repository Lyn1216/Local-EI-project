import json
from motify_functions import *

file_paths1 = ["cellcollective/Arabidopsis thaliana Cell Cycle_26340681.txt",
              "cellcollective/Trichostrongylus retortaeformis_22253585.txt",
              "cellcollective/B bronchiseptica and T retortaeformis coinfection_22253585.txt",
              "cellcollective/CD4+ T cell Differentiation_22871178.txt",
              "cellcollective/Death Receptor Signaling_20221256.txt",
              "cellcollective/Colitis-associated colon cancer_26446703.txt",
              "cellcollective/Signaling in Macrophage Activation_18433497.txt",
              "cellcollective/T cell differentiation_16542429.txt"]
file_paths2 = ["cellcollective/Tumour Cell Invasion and Migration_26528548.txt",
              "cellcollective/Apoptosis Network_19422837.txt",
              "cellcollective/Aurora Kinase A in Neuroblastoma_26616283.txt",
              "cellcollective/Body Segmentation in Drosophila 2013_23520449.txt",
              "cellcollective/B bronchiseptica and T retortaeformis coinfection_22253585.txt",
              "cellcollective/Bordetella bronchiseptica_22253585.txt"]
file_paths3 = ["cellcollective/Apoptosis Network_19422837.txt"]
print("*"*100)
print("FBL")
for file_path in file_paths1:
    graph = create_digraph_from_bn(file_path)
    subgraphs, _ = FBL(file_path, graph, 3)
    name_to_index = name_to_order(file_path)
    filtered_subgraphs = filter_subgraphs(subgraphs, name_to_index)
    dicts = {}
    for filtered_subgraph in filtered_subgraphs:
        dicts[filtered_subgraph] = [(name_to_index[n1], name_to_index[n2], name_to_index[n3]) 
                                    for n1, n2, n3 in [filtered_subgraph]]
    syns = {}
    for key, value in dicts.items():
        un,un_en,syn=text_bn_graph(textfile=str(file_path), candidate_sys=list(value[0]), 
                                   fill_onenode=True, noise=0, save_onenote=False)
        syns[(key, value[0])] = [un,un_en,syn]
    print(file_path)
    print(syns)

print("*"*100)
print("fanin")
for file_path in file_paths2:
    graph = create_digraph_from_bn(file_path)
    subgraphs, _ = fanin(file_path, graph, 3)
    name_to_index = name_to_order(file_path)
    filtered_subgraphs = filter_subgraphs(subgraphs, name_to_index)
    dicts = {}
    for filtered_subgraph in filtered_subgraphs:
        dicts[filtered_subgraph] = [(name_to_index[n1], name_to_index[n2], name_to_index[n3]) 
                                    for n1, n2, n3 in [filtered_subgraph]]
    syns = {}
    for key, value in dicts.items():
        un,un_en,syn=text_bn_graph(textfile=str(file_path), candidate_sys=list(value[0]), 
                                   fill_onenode=True, noise=0, save_onenote=False)
        syns[(key, value[0])] = [un,un_en,syn]
    print(file_path)
    print(syns)

print("*"*100)
print("mutualin")
for file_path in file_paths3:
    graph = create_digraph_from_bn(file_path)
    subgraphs, _ = mutualin(file_path, graph, 3)
    name_to_index = name_to_order(file_path)
    filtered_subgraphs = filter_subgraphs(subgraphs, name_to_index)
    dicts = {}
    for filtered_subgraph in filtered_subgraphs:
        dicts[filtered_subgraph] = [(name_to_index[n1], name_to_index[n2], name_to_index[n3]) 
                                    for n1, n2, n3 in [filtered_subgraph]]
    syns = {}
    for key, value in dicts.items():
        un,un_en,syn=text_bn_graph(textfile=str(file_path), candidate_sys=list(value[0]), 
                                   fill_onenode=True, noise=0, save_onenote=False)
        syns[(key, value[0])] = [un,un_en,syn]
    print(file_path)
    print(syns)