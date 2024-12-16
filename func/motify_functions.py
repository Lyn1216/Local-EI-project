import sys
sys.path.append("../..")
import func.load_database13 as db
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from func.grn_tpm import text_bn_graph
from pathlib import Path
import itertools
import contextlib
import io
import os

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

# 三节点基序
def fanout(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[1], nodes[0]) and \
        nhe(G, nodes[2], nodes[0]) and nhe(G, nodes[2], nodes[1]) and nhe(G, nodes[1], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-fanout{i+1}" for i in range(len(subgraphs))]
    #print(name_subgraphs)
    return subgraphs, name_subgraphs

def fanin(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and nhe(G, nodes[0], nodes[1]) and \
        nhe(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[1]) and nhe(G, nodes[1], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-fanin{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def cascade(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[2], nodes[0]) and nhe(G, nodes[1], nodes[0]) and \
        nhe(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[1]) and nhe(G, nodes[1], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-cascade{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def mutualout(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[2], nodes[0]) and he(G, nodes[0], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and nhe(G, nodes[2], nodes[1]) and nhe(G, nodes[1], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-mutualout{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def mutualin(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and nhe(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and nhe(G, nodes[2], nodes[0]) and he(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-mutualin{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def bimutual(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[1], nodes[2]) and \
        he(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and nhe(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-bimutual{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def FFL(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and nhe(G, nodes[2], nodes[0]) and nhe(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-FFL{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def FBL(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and nhe(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and nhe(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-FBL{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def regulatingmutual(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and nhe(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        he(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and he(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-regulatingmutual{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def regulatedmutual(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and nhe(G, nodes[2], nodes[0]) and he(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-regulatedmutual{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def mutualcascade(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and nhe(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        nhe(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and he(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-mutualcascade{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def semiclique(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        he(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and nhe(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-semiclique{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def clique(file_path, G, k=3):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[0], nodes[2]) and he(G, nodes[1], nodes[2]) and \
        he(G, nodes[1], nodes[0]) and he(G, nodes[2], nodes[0]) and he(G, nodes[2], nodes[1]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-clique{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

# 四节点基序
def bifan(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and nhe(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and nhe(G, nodes[3], nodes[0]) and he(G, nodes[1], nodes[2]) and nhe(G, nodes[2], nodes[1]) and\
        he(G, nodes[1], nodes[3]) and nhe(G, nodes[3], nodes[1]) and nhe(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-bifan{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def twoFFL(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and nhe(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and nhe(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and nhe(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-twoFFL{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def fourone(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and nhe(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and nhe(G, nodes[3], nodes[1]) and nhe(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-fourone{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def FBL4(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        nhe(G, nodes[0], nodes[3]) and nhe(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and nhe(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-fourtwo{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def fiveone(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and nhe(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and nhe(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and nhe(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-fiveone{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def fivetwo(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and nhe(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and nhe(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-fivetwo{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def fivethree(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and nhe(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and nhe(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-fivethree{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def sixone(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and nhe(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and nhe(G, nodes[2], nodes[1]) and\
        he(G, nodes[1], nodes[3]) and nhe(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-sixone{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def sixtwo(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and nhe(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-sixtwo{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def sixthree(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and nhe(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-sixthree{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def sixfour(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and nhe(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-sixfour{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def sevenone(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and nhe(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-sevenone{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def seventwo(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and nhe(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-seventwo{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def seventhree(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-seventhree{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def eightone(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and he(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and nhe(G, nodes[3], nodes[0]) and nhe(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-eightone{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def eighttwo(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if nhe(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and nhe(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and he(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-eighttwo{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def nine(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and he(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and he(G, nodes[1], nodes[2]) and nhe(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-nine{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def ten(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and he(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and he(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and nhe(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-ten{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def eleven(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and he(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and he(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        nhe(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and he(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-eleven{i+1}" for i in range(len(subgraphs))]
    return subgraphs, name_subgraphs

def twelve(file_path, G, k=4):
    node_permutations = list(itertools.permutations(G.nodes(), k))
    subgraphs = []
    for nodes in node_permutations:
        if he(G, nodes[0], nodes[1]) and he(G, nodes[1], nodes[0]) and he(G, nodes[0], nodes[2]) and he(G, nodes[2], nodes[0]) and\
        he(G, nodes[0], nodes[3]) and he(G, nodes[3], nodes[0]) and he(G, nodes[1], nodes[2]) and he(G, nodes[2], nodes[1]) and\
        he(G, nodes[1], nodes[3]) and he(G, nodes[3], nodes[1]) and he(G, nodes[2], nodes[3]) and he(G, nodes[3], nodes[2]):
            subgraphs.append(nodes)
    subgraphs = list(set(tuple(sorted(t)) for t in subgraphs))
    name_subgraphs = [f"{str(file_path)[15:-4]}-twelve{i+1}" for i in range(len(subgraphs))]
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

# def filter_subgraph_dict(subgraphs, name_subgraphs, filtered_subgraphs):
#     subgraph_dict = {name: subgraph for name, subgraph in zip(name_subgraphs, subgraphs)}
#     filtered_subgraph_dict = {name: subgraph for name, subgraph in subgraph_dict.items() if subgraph in filtered_subgraphs}
#     return filtered_subgraph_dict

# 计算un,un_en，syn的值
def cal_motify(function, file_path, n):
    # 文件-图
    graph = create_digraph_from_bn(file_path)
    # 找出三元完全子图
    subgraphs, name_subgraphs = function(file_path, graph, n)
    # 将名称改为序号
    name_to_index = name_to_order(file_path)
    # 排除包含环境的三元组
    filtered_subgraphs = filter_subgraphs(subgraphs, name_to_index)
    dis = {}
    # 将找到的子图序号元组名称换成序号
    for filtered_subgraph in filtered_subgraphs:
        dis[filtered_subgraph] = [(name_to_index[n1], name_to_index[n2], name_to_index[n3], name_to_index[n4]) 
                                  for n1, n2, n3, n4 in [filtered_subgraph]]
    # 存储结果
    dicts = {}
    un = {}
    un_en = {}
    syn = {}
    for (key, value), name_subgraph in zip(dis.items(), name_subgraphs):
    # for key, value in dicts.items():
        unv,un_env,synv=text_bn_graph(textfile=str(file_path), candidate_sys=list(value[0]), 
                                   fill_onenode=True, noise=0, save_onenote=False)
        un[name_subgraph] = unv
        un_en[name_subgraph] = un_env
        syn[name_subgraph] = synv
        dicts[(key, value[0])] = [unv,un_env,synv]
    # # 将找到的三元组名称换成序号
    # replaced_subgraphs = [(name_to_index[n1], name_to_index[n2], name_to_index[n3]) for n1, n2, n3 in filtered_subgraphs]
    # # 存储结果
    # un = {}
    # un_en = {}
    # syn = {}
    # for replaced_subgraph,name_subgraph in zip(replaced_subgraphs,name_subgraphs):
    #     un[name_subgraph],un_en[name_subgraph],syn[name_subgraph]=\
    #     text_bn_graph(textfile=str(file_path), candidate_sys=list(replaced_subgraph), fill_onenode=True, noise=0, save_onenote=False)
    return un, un_en, syn, dicts

# 绘图函数
def make_one_figure(x_values, y_un, y_un_en, y_syn, figurename):
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

def make_un_unen_figure(y_un, y_un_en, figurename):
    # 绘制un和un_en二维散点图
    plt.figure(figsize=(10, 6)) 
    plt.scatter(y_un, y_un_en, marker='o')
    plt.xlabel('un')
    plt.ylabel('un_en')
    plt.title(f"{str(figurename)}_2un_un_en")
    plt.savefig(f"{str(figurename)}_2un_un_en.png")
    #plt.show()

def make_un_syn_figure(y_un, y_syn, figurename):
    # 绘制un和syn二维散点图
    plt.figure(figsize=(10, 6)) 
    plt.scatter(y_un, y_syn, marker='o')
    plt.xlabel('un')
    plt.ylabel('syn')
    plt.title(f"{str(figurename)}_2un_syn")
    plt.savefig(f"{str(figurename)}_2un_syn.png")
    #plt.show()

def make_un_en_syn_figure(y_un_en, y_syn, figurename):
    # 绘制syn和un_en二维散点图
    plt.figure(figsize=(10, 6)) 
    plt.scatter(y_syn, y_un_en, marker='o')
    plt.xlabel('syn')
    plt.ylabel('un_en')
    plt.title(f"{str(figurename)}_2syn_un_en")
    plt.savefig(f"{str(figurename)}_2syn_un_en.png")
    #plt.show()

def make_3D_figure(y_un, y_un_en, y_syn, figurename):
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

def make_un_distribution_figure(y_un, figurename):
    # 绘制un直方图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(y_un, bins=10, color='blue', alpha=0.7, label='UN')
    ax.set_xlabel('UN Values')
    ax.set_ylabel('Frequency')
    ax.set_title(f"{str(figurename)}_UN Distribution")
    plt.savefig(f"{str(figurename)}_UN Distribution.png")
    #plt.show()

def make_un_en_distribution_figure(y_un_en, figurename):
    # 绘制un_en直方图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(y_un_en, bins=10, color='blue', alpha=0.7, label='UN_EN')
    ax.set_xlabel('UN_EN Values')
    ax.set_ylabel('Frequency')
    ax.set_title(f"{str(figurename)}_UN_EN Distribution")
    plt.savefig(f"{str(figurename)}_UN_EN Distribution.png")
    #plt.show()

def make_syn_distribution_figure(y_syn, figurename):
    # 绘制syn直方图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(y_syn, bins=10, color='blue', alpha=0.7, label='SYN')
    ax.set_xlabel('SYN Values')
    ax.set_ylabel('Frequency')
    ax.set_title(f"{str(figurename)}_SYN Distribution")
    plt.savefig(f"{str(figurename)}_SYN Distribution.png")
    #plt.show()

def make_un_un_en_distribution_figure(y_un, y_un_en, figurename):
    # 绘制un和un_en直方图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist([y_un, y_un_en], bins=10, color=["blue", "green"], alpha=1, label=['un',"un_en"], align='mid')
    ax.set_xlabel('Values')
    ax.set_ylabel('Frequency')
    ax.set_title(f"{str(figurename)}_UN&UN_EN Distribution")
    plt.savefig(f"{str(figurename)}_UN&UN_EN Distribution.png")
    ax.legend() 
    #plt.show()

def make_un_bar_figure(un, figurename, motifyname):
   # motifyname = ["fanout", "fanin", "cascade", "mutualout", "mutualin", "bimutual", 
    #              "FFL", "FBL", "regulatingmutual", "regulatedmutual", "mutualcascade", "semiclique", "clique"]
    grouped_values = {motif: [] for motif in motifyname}
    for key, value in un.items():
        for motif in motifyname:
            if motif in key:
                grouped_values[motif].append(value)
                break
    means = [np.mean(grouped_values[motif]) for motif in motifyname]
    stds = [np.std(grouped_values[motif]) for motif in motifyname]

    # 绘制柱状图并添加标准差
    x = np.arange(len(motifyname))
    width = 0.35

    fig, ax = plt.subplots(figsize=(15, 10))
    bars = ax.bar(x, means, width, yerr=stds, capsize=5, color='#74B38F')

    # 添加标题和标签
    #ax.set_title(f"dun about {str(figurename)}")
    ax.set_xlabel('Motifs')
    ax.set_ylabel('Un_sys Values')
    ax.set_xticks(x)
    ax.set_xticklabels(motifyname)
    plt.xticks(rotation=30)
    # ax.legend()
    plt.savefig(f"{str(figurename)}_un_sys Mean and Standard Deviation.png", dpi=600)
    # 显示图形
    #plt.show()

def make_un_en_bar_figure(un_en, figurename, motifyname):
    #motifyname = ["fanout", "fanin", "cascade", "mutualout", "mutualin", "bimutual", 
    #         "FFL", "FBL", "regulatingmutual", "regulatedmutual", "mutualcascade", "semiclique", "clique"]
    grouped_values = {motif: [] for motif in motifyname}
    for key, value in un_en.items():
        for motif in motifyname:
            if motif in key:
                grouped_values[motif].append(value)
                break
    means = [np.mean(grouped_values[motif]) for motif in motifyname]
    stds = [np.std(grouped_values[motif]) for motif in motifyname]

    # 绘制柱状图并添加标准差
    x = np.arange(len(motifyname))
    width = 0.35

    fig, ax = plt.subplots(figsize=(15, 10))
    bars = ax.bar(x, means, width, yerr=stds, capsize=5, color='#9B76B2')

    # 添加标题和标签
    #ax.set_title(f"dun_en about {str(figurename)}")
    ax.set_xlabel('Motifs')
    ax.set_ylabel('Un_en Values')
    ax.set_xticks(x)
    ax.set_xticklabels(motifyname)
    plt.xticks(rotation=30)
    #ax.legend()
    plt.savefig(f"{str(figurename)}_un_en Mean and Standard Deviation.png", dpi=600)
    # 显示图形
    #plt.show()

def make_syn_bar_figure(syn, figurename, motifyname):
    #motifyname = ["fanout", "fanin", "cascade", "mutualout", "mutualin", "bimutual", 
    #         "FFL", "FBL", "regulatingmutual", "regulatedmutual", "mutualcascade", "semiclique", "clique"]
    grouped_values = {motif: [] for motif in motifyname}
    for key, value in syn.items():
        for motif in motifyname:
            if motif in key:
                grouped_values[motif].append(value)
                break
    means = [np.mean(grouped_values[motif]) for motif in motifyname]
    stds = [np.std(grouped_values[motif]) for motif in motifyname]

    # 绘制柱状图并添加标准差
    x = np.arange(len(motifyname))
    width = 0.35

    fig, ax = plt.subplots(figsize=(15, 10))
    bars = ax.bar(x, means, width, yerr=stds, capsize=5, color='#2A69B3')

    # 添加标题和标签
    #ax.set_title(f"dsyn about {str(figurename)}")
    ax.set_xlabel('Motifs')
    ax.set_ylabel('Syn Values')
    ax.set_xticks(x)
    ax.set_xticklabels(motifyname)
    plt.xticks(rotation=30)
    #ax.legend()
    plt.savefig(f"{str(figurename)}_syn Mean and Standard Deviation.png", dpi=600)
    # 显示图形
    #plt.show()

def make_combined_bar_figure(un, un_en, syn, figurename, motifyname):
    #motifyname = ["fanout", "fanin", "cascade", "mutualout", "mutualin", "bimutual", 
    #         "FFL", "FBL", "regulatingmutual", "regulatedmutual", "mutualcascade", "semiclique", "clique"]
    
    # 分组数据并计算均值和标准差
    def group_data(data):
        grouped_values = {motif: [] for motif in motifyname}
        for key, value in data.items():
            for motif in motifyname:
                if motif in key:
                    grouped_values[motif].append(value)
        means = [np.mean(grouped_values[motif]) for motif in motifyname]
        stds = [np.std(grouped_values[motif]) for motif in motifyname]
        return means, stds

    un_means, un_stds = group_data(un)
    un_en_means, un_en_stds = group_data(un_en)
    syn_means, syn_stds = group_data(syn)

    x = np.arange(len(motifyname))
    width = 0.25  # 减小宽度以便三个柱子可以并排显示

    # 创建图表和第一个轴
    fig, ax1 = plt.subplots(figsize=(15, 10))

    # 绘制syn数据集的柱状图
    bars_syn = ax1.bar(x - width, syn_means, width, yerr=syn_stds, error_kw={'capsize': 3, 'elinewidth': 1, 'ecolor': 'grey'}, color='#2A69B3', label='Syn Values')
    ax1.set_ylim(0,0.7) 
    # 创建第二个轴对象
    ax2 = ax1.twinx()

    # 绘制un和un_en数据集的柱状图
    bars_un = ax2.bar(x, un_means, width, yerr=un_stds, error_kw={'capsize': 3, 'elinewidth': 1, 'ecolor': 'grey'}, color='#74B38F', label='Un_sys Values', alpha=0.7)
    bars_un_en = ax2.bar(x + width, un_en_means, width, yerr=un_en_stds, error_kw={'capsize': 3, 'elinewidth': 1, 'ecolor': 'grey'}, color='#9B76B2', label='Un_en Values', alpha=0.7)
    ax2.set_ylim(0,2.5) 
    # 设置标题和标签
    ax1.set_xlabel('Motifs')
    ax1.set_ylabel('Syn Values')
    ax2.set_ylabel('Un_sys and Un_en Values')

    ax1.set_xticks(x)
    ax1.set_xticklabels(motifyname, rotation=30)

    # 设置图例
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc='upper right')

    # 保存图表
    plt.savefig(f"{str(figurename)}_combined Mean and Standard Deviation.png", dpi=600)
