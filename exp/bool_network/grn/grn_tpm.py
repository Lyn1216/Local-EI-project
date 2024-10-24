import load_database13 as db
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from itertools import product

def tpm_one(f_one, inputs, ss, noise):
    lens = len(f_one)
    matrix = np.zeros([lens, 2])
    for i in range(lens):
        if f_one[i] == 0:
            matrix[i, 0] = 1 - noise
        else:
            matrix[i, 0] = noise
        matrix[i, 1] = 1 - matrix[i, 0]
        
    if ss not in inputs:
        matrix = np.tile(matrix, (2, 1))
        en_size = len(inputs)
    else:
        en_size = len(inputs) - 1
    return matrix, en_size
    
def text_bn_graph(textfile = 'example.txt', condidate_sys=None, fill_onenode=False, noise=0):
    F, I, degree, variables, constants = db.text_to_BN(folder='',textfile=textfile)
    
    G = nx.DiGraph()
    all_nodes = variables + constants
    G.add_nodes_from(all_nodes)
    
    for i in range(len(variables)):
        for j in I[i]:
            G.add_edges_from([(all_nodes[j], variables[i])])
    
    pos = nx.spring_layout(G)  # 为图形设置布局
    nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=700, edge_color='k', linewidths=1, font_size=15)
    plt.show()
    print("all intrinsic variables: " + ','.join(variables))
    print("external parameters:" + ','.join(constants))
    if fill_onenode is False:
        for i in range(len(variables)):
            tpm1, en_size = tpm_one(F[i], I[i], i, noise=noise)
            print("mechanism:    " + variables[i])
            print("environment:    " + ','.join([all_nodes[j] for j in I[i]]))
            print(tpm1)
            print("un:  " + str(unique(tpm1, 1, en_size)[0]))
            print("syn:  " + str(synergy(tpm1, 1, en_size)))
            print("vividness:  " + str(condi_ei(tpm1, 1, en_size)))
            print(120 * '-')

    if condidate_sys is not None:
        neigbors, tpm = tpm_comb(condidate_sys, F, I, noise)
        print("mechanism:    " + ','.join([variables[j] for j in condidate_sys]))
        print("environment:    " + ','.join([all_nodes[j] for j in neigbors]))
        print(tpm)
        print("un:  " + str(unique(tpm, len(condidate_sys), len(neigbors)-len(condidate_sys))[0]))
        print("syn:  " + str(synergy(tpm, len(condidate_sys), len(neigbors)-len(condidate_sys))))
        print("vividness:  " + str(condi_ei(tpm, len(condidate_sys), len(neigbors)-len(condidate_sys))))
        print(120 * '-')

def permute_matrix_rows(original_order, new_order):
    # 原始和新顺序的长度（应该是3）
    n = len(original_order)
    # 创建一个映射字典
    mapping = []
    
    # 生成所有可能的状态（000到111）
    for state in product([0, 1], repeat=n):
        # 将元组转换为二进制字符串
        original_state = ''.join(map(str, state))
        #new_state = ''.join(original_state[i] for i in [original_order.index(new_order[j]) for j in range(n)])
        # 根据新顺序重新排列状态位
        #indices = np.where(original_order == new_order[j] for j in range(n))[0][0]  # 注意：这里假设每个元素都是唯一的

        # 然后你可以使用这些索引来重新排列状态位
        new_state = ''.join(original_state[i] for i in [np.where(original_order == new_order[j])[0][0] for j in range(n)])
        
        # 将原始状态映射到新状态
        mapping.append(new_state)
    
    return mapping

def add_missing_elements(array1, array2):
    # 将数组转换为集合
    set1 = set(array1)
    set2 = set(array2)
    array22 = array2.copy()
    # 找出在set1中但不在set2中的元素
    missing_elements = set1 - set2
    
    # 将缺失的元素添加到array2中
    array22 = np.concatenate((np.array(list(missing_elements)), array22))
    
    return array22

def decimal_to_binary(decimal, min_length=1):
    if decimal == 0:
        return "0" if min_length == 1 else "0".zfill(min_length)
    binary = ""
    while decimal > 0:
        binary = str(decimal % 2) + binary
        decimal = decimal // 2
    # 使用 zfill 确保二进制字符串至少有 min_length 长度
    return binary.zfill(min_length)

def tpm_comb(candidate_system, F, I, noise):
    neigbors = candidate_system
    tpm_list = []
    for i in candidate_system:
        neigbors = np.concatenate((neigbors, I[i]))
        tpm_list.append(tpm_one(F[i], I[i], i, noise=noise)[0])
    neigbors = np.unique(neigbors)
    matrix = np.ones([2**len(neigbors), 2**len(candidate_system)])
    for i, s in enumerate(candidate_system):
        new_arr = add_missing_elements(neigbors, I[s])
        while tpm_list[i].shape[0] < matrix.shape[0]:
            tpm_list[i] = np.tile(tpm_list[i], (2, 1))
        mapping = permute_matrix_rows(new_arr, neigbors)
        new_id = [int(index,2) for index in mapping]
        tpm_list[i] = tpm_list[i][new_id]
    for i in range(matrix.shape[1]):
        ind_str = decimal_to_binary(i, min_length=len(candidate_system))
        for k, j in enumerate(ind_str):
            matrix[:, i] *= tpm_list[k][:, int(j)] 
    return neigbors, matrix



def tpm_ei(tpm, log_base = 2):
    # marginal distribution of y given x ~ Unifrom Dist
    puy = tpm.sum(axis=0)
    n = tpm.shape[0]
    # replace 0 to a small positive number to avoid log error
    eps = 1E-10
    tpm_e = np.where(tpm==0, eps, tpm)
    puy_e = np.where(tpm==0, eps, puy)
    
    # calculate EI of specific x
    ei_x = (np.log2(n * tpm_e / puy_e) / np.log2(log_base)  * tpm).sum(axis=1)
    
    # calculate total EI
    ei_all = ei_x.mean()
    return ei_all


def synergy(markov_matrix, mech_size, en_size):
    ei_all = condi_ei(markov_matrix, mech_size, en_size)
    un = unique(markov_matrix, mech_size, en_size)[0]
    syn = ei_all - un
    return syn

def condi_ei(markov_matrix, mech_size, en_size):
    ei = 0
    state_size = 2**mech_size
    state_en = 2**en_size
    for e in range(state_en):
        local_markov = np.zeros([state_size, state_size])
        for num in range(state_size):
            binary_string = decimal_to_binary(num, min_length=mech_size)
            padded_binary_string = binary_string + decimal_to_binary(e, min_length=en_size)
            binary_array = [int(bit) for bit in padded_binary_string] 
            pattern = int(''.join(str(cell) for cell in binary_array), 2)
            local_markov[num, :] = markov_matrix[pattern, :]
        ei += tpm_ei(local_markov)
    ei = ei / state_en
    return ei 

def unique(markov_matrix, mech_size, en_size):
    state_size = 2**mech_size
    state_en = 2**en_size
    mixed_markov = np.zeros([state_size, state_size])
    for e in range(state_en):
        local_markov = np.zeros([state_size, state_size])
        for num in range(state_size):
            binary_string = decimal_to_binary(num, min_length=mech_size)
            padded_binary_string = binary_string + decimal_to_binary(e, min_length=en_size)
            binary_array = [int(bit) for bit in padded_binary_string] 
            pattern = int(''.join(str(cell) for cell in binary_array), 2)
            local_markov[num, :] = markov_matrix[pattern, :]
        mixed_markov += local_markov
    ei = tpm_ei(mixed_markov / state_en)
    return ei, mixed_markov