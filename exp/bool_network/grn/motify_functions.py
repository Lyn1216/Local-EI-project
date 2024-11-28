import sys
sys.path.append("../..")
import func.load_database13 as db
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from grn_tpm import text_bn_graph
from pathlib import Path
import itertools
import contextlib
import io

# 生成图的函数
def create_digraph_from_bn(file_path):
    F, I, degree, variables, constants = db.text_to_BN(folder="", textfile=str(file_path))
    G = nx.DiGraph()
    all_nodes = variables + constants
    G.add_nodes_from(all_nodes)
    for i in range(len(variables)):
        for j in I[i]:
            G.add_edges_from([(all_nodes[j], variables[i])])
    return G

# 两节点之间有边
def he(G, u, v):
    return G.has_edge(u, v)

# 两节点之间无边
def nhe(G, u, v):
    return not G.has_edge(u, v)

def fanout(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[1], nodes[0]) and \
        nhe(G, nodes[2], nodes[0]) and nhe(G, nodes[2], nodes[1]) and nhe(G, nodes[1], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-fanout{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def fanin(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and nhe(G, nodes[0], nodes[1]) and \
        nhe(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[1]) and nhe(G, nodes[1], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-fanin{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def cascade(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[2], nodes[0]) and nhe(G, nodes[1], nodes[0]) and \
        nhe(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[1]) and nhe(G, nodes[1], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-cascade{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def mutualout(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[2], nodes[0]) and he(G, nodes[0], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and nhe(G, nodes[2], nodes[1]) and nhe(G, nodes[1], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-mutualout{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def mutualin(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and nhe(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and nhe(G, nodes[2], nodes[0]) and he(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-mutualin{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def bimutual(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[1], nodes[2]) and \
        he(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and nhe(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-bimutual{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def FFL(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and nhe(G, nodes[2], nodes[0]) and nhe(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-FFL{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def FBL(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and nhe(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and nhe(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-FBL{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def regulatingmutual(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and nhe(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        he(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and he(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-regulatingmutual{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def regulatedmutual(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and nhe(G, nodes[2], nodes[0]) and he(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-regulatedmutual{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def mutualcascade(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and nhe(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and he(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-mutualcascade{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def semiclique(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        he(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and nhe(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-semiclique{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def clique(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        he(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and he(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[:6]}-clique{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

# 将文件的节点名称与序号对应
def name_to_order(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    name_to_index = {}
    for index, line in enumerate(lines):
        line = line.strip()
        equal_sign_index = line.find('=')
        if equal_sign_index != -1:
            content = line[:equal_sign_index].strip()
            name_to_index[content] = index
    return name_to_index

# 去除包含环境的子图
def filter_subgraphs(subgraphs, name_to_index):
    # 使用列表推导式过滤子图
    filtered_subgraphs = [
        subgraph for subgraph in subgraphs
        if all(name in name_to_index for name in subgraph)
    ]
    return filtered_subgraphs

# 计算un,un_en，syn的值
def cal_motify(function, file_path):
    # 文件-图
    graph = create_digraph_from_bn(file_path)
    # 找出三元完全子图
    subgraphs, name_subgraphs = function(file_path, graph, 3)
    # 将名称改为序号
    name_to_index = name_to_order(file_path)
    # 排除包含环境的三元组
    filtered_subgraphs = filter_subgraphs(subgraphs, name_to_index)
    # 将找到的三元组名称换成序号
    replaced_subgraphs = [(name_to_index[n1], name_to_index[n2], name_to_index[n3]) for n1, n2, n3 in filtered_subgraphs]
    # 存储结果
    un = {}
    un_en = {}
    syn = {}
    for replaced_subgraph,name_subgraph in zip(replaced_subgraphs,name_subgraphs):
        un[name_subgraph],un_en[name_subgraph],syn[name_subgraph]=\
        text_bn_graph(textfile=str(file_path), candidate_sys=list(replaced_subgraph), fill_onenode=True, noise=0, save_onenote=False)
    return un, un_en, syn

# 绘图函数
def make_figure(x_values, y_un, y_un_en, y_syn, figurename):
    # 绘制一维散点图
    plt.figure(figsize=(10, 6)) 
    plt.scatter(x_values, y_un, label="un", marker='o')
    plt.scatter(x_values, y_un_en, label="un_en", marker='o')
    plt.scatter(x_values, y_syn, label="syn", marker='o')
    plt.legend()   
    plt.xlabel('Subgraphs')
    plt.ylabel('Values')
    plt.title(f"{str(figurename)}_one")
    plt.savefig(f"{str(figurename)}_one.png")
    #plt.show()
    # 绘制un和un_en二维散点图
    plt.figure(figsize=(10, 6)) 
    plt.scatter(y_un, y_un_en, marker='o')
    plt.xlabel('un')
    plt.ylabel('un_en')
    plt.title(f"{str(figurename)}_2un_un_en")
    plt.savefig(f"{str(figurename)}_2un_un_en.png")
    #plt.show()
    # 绘制un和syn二维散点图
    plt.figure(figsize=(10, 6)) 
    plt.scatter(y_un, y_syn, marker='o')
    plt.xlabel('un')
    plt.ylabel('syn')
    plt.title(f"{str(figurename)}_2un_syn")
    plt.savefig(f"{str(figurename)}_2un_syn.png")
    #plt.show()
    # 绘制syn和un_en二维散点图
    plt.figure(figsize=(10, 6)) 
    plt.scatter(y_syn, y_un_en, marker='o')
    plt.xlabel('syn')
    plt.ylabel('un_en')
    plt.title(f"{str(figurename)}_2syn_un_en")
    plt.savefig(f"{str(figurename)}_2syn_un_en.png")
    #plt.show()
    # 绘制三维散点图
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(y_un, y_un_en, y_syn, marker='o')
    ax.set_xlabel('UN')
    ax.set_ylabel('UN_EN')
    ax.set_zlabel('SYN')
    ax.set_title(f"{str(figurename)}_3D")
    plt.savefig(f"{str(figurename)}_3D.png")
    #plt.show()
    # 绘制un直方图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(y_un, bins=10, color='blue', alpha=0.7, label='UN')
    ax.set_xlabel('UN Values')
    ax.set_ylabel('Frequency')
    ax.set_title(f"{str(figurename)}_UN Distribution")
    plt.savefig(f"{str(figurename)}_UN Distribution.png")
    #plt.show()
    # 绘制un_en直方图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(y_un_en, bins=10, color='blue', alpha=0.7, label='UN_EN')
    ax.set_xlabel('UN_EN Values')
    ax.set_ylabel('Frequency')
    ax.set_title(f"{str(figurename)}_UN_EN Distribution")
    plt.savefig(f"{str(figurename)}_UN_EN Distribution.png")
    #plt.show()
    # 绘制syn直方图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(y_syn, bins=10, color='blue', alpha=0.7, label='SYN')
    ax.set_xlabel('SYN Values')
    ax.set_ylabel('Frequency')
    ax.set_title(f"{str(figurename)}_SYN Distribution")
    plt.savefig(f"{str(figurename)}_SYN Distribution.png")
    #plt.show()
    # 绘制un和un_en直方图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist([y_un, y_un_en], bins=10, color=["blue", "green"], alpha=1, label=['un',"un_en"], align='mid')
    ax.set_xlabel('Values')
    ax.set_ylabel('Frequency')
    ax.set_title(f"{str(figurename)}_UN&UN_EN Distribution")
    plt.savefig(f"{str(figurename)}_UN&UN_EN Distribution.png")
    ax.legend() 
    #plt.show()

