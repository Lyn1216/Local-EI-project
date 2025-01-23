# import sys
# sys.path.append("../../..")
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from itertools import product
from func.EI_calculation import tpm_ei_new, tpm_ei_new2
import func.load_database13 as db

def tpm_one(f_one, inputs, ss, noise):
#     if f_one == []:
#         print("There is an empty F")
#         f_one = [0, 1]
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
        all_inputs = np.concatenate(([ss], inputs))
    else:
        en_size = len(inputs) - 1
        all_inputs = inputs
    return matrix, en_size, all_inputs

def text_bn_graph(folder = '', textfile = 'example.txt', candidate_sys=None, figure_show=False, noise=0):
    F, I, degree, variables, constants = db.text_to_BN(folder=folder,textfile=textfile)
    if candidate_sys == "all":
        candidate_sys = range(len(variables))

    all_nodes = variables + constants
    
    if figure_show:
        G = nx.DiGraph()
        G.add_nodes_from(all_nodes)
        for i in range(len(variables)):
            for j in I[i]:
                G.add_edges_from([(all_nodes[j], variables[i])])

        pos = nx.spring_layout(G)  # 为图形设置布局
        nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=700, edge_color='k', linewidths=1, font_size=15)
        plt.show()
        print("all intrinsic variables: " + ','.join(variables))
        print("external parameters:" + ','.join(constants))
#     onenote_tpm_result = {}
#     onenote_un_result = {}
#     onenote_syn_result = {}
#     onenote_vividness_result = {}
#     if save_onenote is True :
#         for i in range(len(variables)):
#             tpm1, en_size, _ = tpm_one(F[i], I[i], i, noise=noise)
#             onenote_tpm_result[variables[i]] = tpm1
#             onenote_un_result[variables[i]] = unique(tpm1, 1, en_size)[0]
#             onenote_syn_result[variables[i]] = synergy(tpm1, 1, en_size)[0]
#             onenote_vividness_result[variables[i]] = condi_ei(tpm1, 1, en_size)
#             if fill_onenode  is False:
#                 print("mechanism:    " + variables[i])
#                 print("environment:    " + ','.join([all_nodes[j] for j in I[i]]))
#                 print(tpm1)
#                 print("un:  " + str(unique(tpm1, 1, en_size)[0]))
#                 print("syn:  " + str(synergy(tpm1, 1, en_size)[0]))
#                 #print("vividness:  " + str(condi_ei(tpm1, 1, en_size)))
#                 print(120 * '-')
#             else:
#                 continue

    if candidate_sys is not None:
#         print("mechanism:    " + ','.join([variables[j] for j in candidate_sys]))
        if len(candidate_sys) <= 20:
            neigbors, tpm = tpm_comb(candidate_sys, F, I, noise)
#             print("tpm: ")
#             print(tpm)
#             print("environment:    " + ','.join([all_nodes[j] for j in neigbors]))
#             un = unique(tpm, len(candidate_sys), len(neigbors)-len(candidate_sys))[0]
#             un_en = en_unique(tpm, len(candidate_sys), len(neigbors)-len(candidate_sys))[0]
            syn, expansive, introverted, tpm_dic = synergy(tpm, len(candidate_sys), len(neigbors)-len(candidate_sys))
        if figure_show:
            print("mechanism:    " + ','.join([variables[j] for j in candidate_sys]))
#             print("un_en:  " + str(un_en))
#             un_approx = un_comb(candidate_sys, F, I, noise)
#             syn_approx = syn_comb(candidate_sys, F, I, noise)
#             print("un_approx:  " + str(un_approx))
#             print("syn_approx:  " + str(syn_approx))
        
#         else:
#             neigbors = nei_comb(candidate_sys, F, I, noise)[0]
#             tpm = "None"
# #             print("environment:    " + ','.join([all_nodes[j] for j in neigbors]))
#             un = un_comb(candidate_sys, F, I, noise)
#             syn = syn_comb(candidate_sys, F, I, noise)

#         print("un:  " + str(un))
#         print("syn:  " + str(syn))
#         #print("vividness:  " + str(vivid))
#         #condi_ei(tpm, len(condidate_sys), len(neigbors)-len(condidate_sys))
#         print(120 * '-')

        return expansive, introverted, syn, tpm, len(neigbors)



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
    if min_length == 0:
        return ''
    if decimal == 0:
        return "0" if min_length == 1 else "0".zfill(min_length)
    binary = ""
    while decimal > 0:
        binary = str(decimal % 2) + binary
        decimal = decimal // 2
    # 使用 zfill 确保二进制字符串至少有 min_length 长度
    return binary.zfill(min_length)

def nei_comb(candidate_system, F, I, noise):
    neigbors = candidate_system
    tpm_list = []
    all_in_list = []
    for i in candidate_system:
        neigbors = np.concatenate((neigbors, I[i]))
        tpm, size, ins = tpm_one(F[i], I[i], i, noise=noise)
        tpm_list.append(tpm)
        all_in_list.append(ins)
    #neigbors = np.unique(neigbors)
    seen = set()
    neigbors_un = [x for x in neigbors if not (x in seen or seen.add(x))]
#     print(neigbors_un)
    neigbors_ = np.setdiff1d(neigbors, candidate_system)
    return neigbors_un, tpm_list, neigbors_, all_in_list

def tpm_comb(candidate_system, F, I, noise):
    neigbors, tpm_list, _, _ = nei_comb(candidate_system, F, I, noise)
    matrix = np.ones([2**len(neigbors), 2**len(candidate_system)])
    for i, s in enumerate(candidate_system):
        new_arr = add_missing_elements(neigbors, I[s])
#         print("new_arr")
#         print(new_arr)
        while tpm_list[i].shape[0] < matrix.shape[0]:
            tpm_list[i] = np.tile(tpm_list[i], (2, 1))
        mapping = permute_matrix_rows(neigbors, new_arr)
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


def synergy(markov_matrix, mech_size, en_size, new_item=False):
    eis, dets, nondegs, tpm_dic  = condi_ei(markov_matrix, mech_size, en_size)
    un, det, nondeg, mixed_markov = unique(markov_matrix, mech_size, en_size)
    syn = eis - un
    if new_item:
        syn_nondeg = nondegs - det - nondeg
        syn_det = dets
        return syn, syn_nondeg, syn_det, tpm_dic
    else:
        expansive = nondegs - det
        introverted = dets - nondeg
        return syn, expansive, introverted, tpm_dic

def condi_ei(markov_matrix, mech_size, en_size):
    eis = 0
    dets = 0
    nondegs = 0
    state_size = 2**mech_size
    state_en = 2**en_size
    tpm_dic = {}
    for e in range(state_en):
        local_markov = np.zeros([state_size, state_size])
        en_str = decimal_to_binary(e, min_length=en_size)
        for num in range(state_size):
            binary_string = decimal_to_binary(num, min_length=mech_size)
            padded_binary_string = binary_string + en_str
            binary_array = [int(bit) for bit in padded_binary_string] 
            pattern = int(''.join(str(cell) for cell in binary_array), 2)
            local_markov[num, :] = markov_matrix[pattern, :]
        ei,det,nondeg = tpm_ei_new2(local_markov)
        eis += ei
        dets += det
        nondegs += nondeg
        tpm_dic[en_str] = local_markov
    eis = eis / state_en
    dets /= state_en
    nondegs /= state_en
    return eis, dets, nondegs, tpm_dic 

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
    ei,det,nondeg = tpm_ei_new2(mixed_markov / state_en)
    return ei, det, nondeg, mixed_markov

def en_unique(markov_matrix, mech_size, en_size):
    state_size = 2**mech_size
    state_en = 2**en_size
    mixed_markov = np.zeros([state_en, state_size])
    for s in range(state_size):
        local_markov = np.zeros([state_en, state_size])
        for num in range(state_en):
            binary_string = decimal_to_binary(num, min_length=en_size)
            padded_binary_string = decimal_to_binary(s, min_length=state_size) + binary_string
            binary_array = [int(bit) for bit in padded_binary_string] 
            pattern = int(''.join(str(cell) for cell in binary_array), 2)
            local_markov[num, :] = markov_matrix[pattern, :]
        mixed_markov += local_markov
    ei,det,deg,eff,det_c,deg_c = tpm_ei_new(mixed_markov / state_size)
    return ei, mixed_markov

def partition(tpm, inputs, eu_list):
    in_rest = np.setdiff1d(inputs, eu_list)
#     print(tpm)
#     print(in_rest)
    #res_ind = np.where(inputs==in_rest)[0]
    #print(res_ind)
    rows = 2**len(in_rest)
    tpm_p = np.zeros([rows, 2])
    for i in range(tpm.shape[0]):
        strs = decimal_to_binary(i, min_length=len(inputs))
        #print(np.where(inputs==in_rest[0]))
        str_res = [strs[np.where(inputs==k)[0][0]] for k in in_rest]
        str_res = ''.join(str_res)
#         print(str_res)
        num = int(str_res, 2)
        tpm_p[num,:] += tpm[i, :]
#         print(tpm_p)
    tpm_p /= 2**len(eu_list)
    return tpm_p

def un_comb(candidate_system, F, I, noise):
    neigbors, tpm_list, neighbors_, all_in_list = nei_comb(candidate_system, F, I, noise)
    un = 0
    set1 = set(neighbors_)
    for i, s in enumerate(candidate_system):
#         print(un)
        set2 = set(I[s])
        intersection = set1.intersection(set2)
        eu_list = list(intersection)
        if eu_list == []:
            un += tpm_ei_new(tpm_list[i])[0]
        else:
            un += tpm_ei_new(partition(tpm_list[i], all_in_list[i], eu_list))[0]
    return un

def syn_comb(candidate_system, F, I, noise):
    neigbors, tpm_list, neighbors_, all_in_list = nei_comb(candidate_system, F, I, noise)
    syn = 0
    set1 = set(neighbors_)
    for i, s in enumerate(candidate_system):
        inputs = all_in_list[i]
        set2 = set(I[s])
        intersection = set1.intersection(set2)
        eu_list = list(intersection)
        if len(eu_list) > 0:
            condi = 0
            in_rest = np.setdiff1d(inputs, eu_list)
            rows = 2**len(in_rest)
            for j in range(2**len(eu_list)):
                tpm_p = np.zeros([rows, 2])
                str_en = decimal_to_binary(j, min_length=len(eu_list))
                for k in range(2**len(inputs)):
                    strs_all = decimal_to_binary(k, min_length=len(inputs))
                    str_eu = [strs_all[np.where(inputs==e)[0][0]] for e in eu_list]
                    str_eu = ''.join(str_eu)
                    if str_eu == str_en:
                        num = int(strs_all, 2)
                        str_in = [strs_all[np.where(inputs==n)[0][0]] for n in in_rest]
                        str_in = ''.join(str_in)
                        num_in = int(str_in, 2)
#                         print("sys:"+str(s))
#                         print("tpm_p:"+str(tpm_p.shape))
#                         print("num_in:"+str(num_in))
#                         print("tpm_all:"+str(tpm_list[i].shape))
#                         print("num:"+str(num))
                        tpm_p[num_in, :] = tpm_list[i][num, :]
                condi += tpm_ei_new(tpm_p)[0]
            syn += condi / 2**len(eu_list)
            syn -= tpm_ei_new(partition(tpm_list[i], all_in_list[i], eu_list))[0]
    return syn

def tpm_to_dis(tpm, mech_size, en_size):
    state_size = 2**mech_size
    state_en = 2**en_size
    tpm_dis = np.zeros([state_size*state_en, state_size])
    for e in range(state_en):
        for num in range(state_size):
            binary_string = decimal_to_binary(num, min_length=mech_size)
            padded_binary_string = binary_string + decimal_to_binary(e, min_length=en_size)
            binary_array = [int(bit) for bit in padded_binary_string] 
            pattern = int(''.join(str(cell) for cell in binary_array), 2)
            tpm_dis[:, num] += tpm[:, pattern]
    return tpm_dis

def iit_tpm_cal(tpm_v, mech_size, en_size, dis=True, new_item=False):
    if dis:
        tpm_dis = tpm_v
    else:
        tpm_dis = tpm_to_dis(tpm_v, mech_size, en_size)
    un = unique(tpm_dis, mech_size, en_size)[0]
    syn, item1, item2, tpm_dic = synergy(tpm_dis, mech_size, en_size, new_item)
    un_en = en_unique(tpm_dis, mech_size, en_size)[0]
#     print("un:  " + str(un))
#     print("un_en:  " + str(un_en))
#     print("syn:  " + str(syn))
    return un, un_en, syn, item1, item2, tpm_dic
    
