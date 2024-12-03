import numpy as np
import pandas as pd
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
from grn_tpm import iit_tpm_cal

def value_boolnet(symbol_boolnet, weights):
    results = {key: weights.get(w, 0) for key, w in symbol_boolnet.items()}
    return results

def get_vars(boolnet):
    """
    Extract variables form Bool net representation
    """
    vars = set()
    for edge in boolnet.keys():
        vars.update(edge)
    assert ''.join(vars).isupper()
    vars = sorted(list(vars))
    return vars


def get_wm(boolnet, vars):
    """
    Get weights matrix
    """
    dim = len(vars)
    # make weights matrix
    wm = pd.DataFrame(np.zeros((dim, dim)), columns=vars, index=vars)
    for s, e in boolnet:
        wm.loc[s, e] = boolnet[(s, e)]
    return wm.values

# def get_wm2(boolnet, vars):
#     """
#     Get weights matrix
#     """
#     dim = len(vars)
#     # make weights matrix
#     wm = np.zeros((dim, dim, dim))
#     for s1, s2, e in boolnet:
#         wm[s1, s2, e] = boolnet[(s1, s2, e)]
#     return wm.values


def get_states(vars, big_endian=True):
    """
    Using upper case and lower case of variable name to represent 2 states of the variable
    and return the full set of combination states of multiple variables
    E.g. ['A', 'B'] -> [('a', 'b'), ('A', 'b'), ('a', 'B'), ('A', 'B') ]
         ['Ab', 'Cd'] -> [('ab', 'cd'), ('AB', 'cd'), ('ab', 'CD'), ('AB', 'CD')]
    @param
    """

    order = 1 if big_endian else -1
    values = [[v.lower(), v.upper()] for v in vars[::order]]
    combinations = list(itertools.product(*values))
    result = [combination[::order] for combination in combinations]
    return result


def state_name(state):
    """
    E.g. ('A', 'B', 'C') -> 'ABC'
    """
    return ''.join(state)


def get_state_values(state):
    """
    assign 1 if all characters of state name are upper case, else -1.
    E.g. ('A', 'B', 'c') -> (1, 1, -1)
         ('AB', 'cd') -> (1, -1)
    """
    values = np.array([1 if v.isupper() else -1 for v in state])
    return values

def make_tpm(bnet, w, k=1, image_show=False, syn_term=False):
    boolnet = value_boolnet(bnet, w)
    vars = list(get_vars(boolnet))
    states = get_states(vars)
    n = len(states)
    s_values = np.array([get_state_values(s) for s in states])
    if syn_term:
        u = s_values.reshape(-1,1)
        uu = u @ u.T
        wm = get_wm2(boolnet, vars)
        states_to_units_pos_tpm = 1. / (1. + np.exp(-k * uu.reshape(-1) @ wm.reshape(n*n,-1)))
    else:
        wm = get_wm(boolnet, vars)
        states_to_units_pos_tpm = 1. / (1. + np.exp(-k * s_values @ wm))
        
    states_to_units_neg_tpm = 1. - states_to_units_pos_tpm
    states_to_units_tpm = np.concatenate([states_to_units_pos_tpm, states_to_units_neg_tpm], axis=1)
    pos_vars = [v.upper() for v in vars]
    neg_vars = [v.lower() for v in vars]
    states_names = [state_name(s) for s in states]
    states_to_units_tpm_df = pd.DataFrame(states_to_units_tpm, columns=pos_vars+neg_vars, index=states_names)
    
    # states: cartesian product of unit probabilitys
    states_tpm = states_to_units_tpm_df.apply(
        lambda row: pd.Series(
            data=[np.array([row[u] for u in s]).prod() for s in states], 
            index=states_names
        ), axis=1)
    if image_show:
        sns.heatmap(states_tpm, annot=True, fmt='.2f', cmap='Greys')
    return states_tpm, states_tpm.values

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

def tpm_series(tpm, init_state, steps, seed = 42):
    np.random.seed(seed)
    init_num = int(init_state, 2)
    serie = [init_num]
    serie_str = [init_state]
    for t in range(steps):
        num = serie[t]
        probabilities = tpm[num, :]
        sample = np.random.choice(range(len(probabilities)), p=probabilities)
        serie.append(sample)
        serie_str.append(decimal_to_binary(sample, min_length=int(np.log2(len(probabilities)))))
    return serie, serie_str

def tpm_syn(k, w):
    w_1b_b = w
    w_2b_b = 0
    w_2c_b = 0
    w_1c_c = 0
    w_2c_c = w
    w_1a_a = 0
    w_2a_a = 0
    w_2b_a = 0
    w_12a_a = w
    w_12b_b = 0
    w_12c_c = 0
    w_a_b = 1 - w
    w_b_c = 1 - w
    w_c_a = 1 - w
    w_a_a = 0
    w_b_b = 0
    w_c_c = 0
    w_c_b = 0
    tpm_list = []
    for j in ['a', 'b', 'c', '1', '2']:
        tpm = np.zeros([32,2])
        for i in range(32):
            inputs = decimal_to_binary(i, min_length=5)
            in_ls = [2*int(n)-1 for n in inputs]
            if j=='a':
                term = w_12a_a*in_ls[0]*in_ls[3]*in_ls[4] + w_2a_a*in_ls[0]*in_ls[4] + w_c_a*in_ls[2] + w_a_a*in_ls[0]
            elif j=='b':
                term = w_1b_b*in_ls[1]*in_ls[3] + w_2b_b*in_ls[1]*in_ls[4] + w_a_b*in_ls[0] + w_b_b*in_ls[1]
            elif j=='c':
                term = w_12c_c*in_ls[2]*in_ls[3]*in_ls[4] + w_2c_c*in_ls[2]*in_ls[4] + w_b_c*in_ls[1] + w_c_c*in_ls[2]
            else:
                term = 0.5

            tpm[i,1] = 1. / (1. + np.exp(-k * term))
            tpm[i,0] = 1 - tpm[i,1]
        tpm_list.append(tpm)

    tpm_all = np.ones([32, 32])
    for i in range(32):
        ind_str = decimal_to_binary(i, min_length=5)
        for m, n in enumerate(ind_str):
            tpm_all[:, i] *= tpm_list[m][:, int(n)] 

    return tpm_all

def serie_plot(bnet, w, k, steps, seed=1, name='', leg=False, syn_inter=False, figure_show=False):
    series_ls = []
    if syn_inter:
        tpm_v = tpm_syn(k=k, w=w)
    else:
        tpm, tpm_v = make_tpm(bnet, w=w, k=k)
    un_sys, un_en, syn, tpm_dic = iit_tpm_cal(tpm_v, mech_size=3, en_size=2)   
    colors = ["#BB4F4F", '#2A69B3', '#74B38F', '#FFA500']
    strs = [decimal_to_binary(i, min_length=3) for i in range(8)]
    for init_state in strs:
        #series_en = []
        if figure_show:
            fig, ax = plt.subplots(figsize=(5,2))
        for indx,en in enumerate(["00", "01", "10", "11"]):
            en_state = en
            serie, serie_str = tpm_series(tpm_dic[en_state], init_state, steps, seed)
            series_ls.append(serie)
            if figure_show:
                ax.scatter(range(steps+1), serie, label='en_state:'+en_state, color=colors[indx])

        if figure_show:
            ax.set_xlabel('Time')
            ax.set_ylabel('System state')

            # 设置y轴的标签
            ax.set_yticks(range(8))
            ax.set_yticklabels(strs)
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            if leg:
                plt.legend(by_label.values(), by_label.keys(), loc=[1.01, 0])
            plt.title(name + '_init=' + init_state + '_syn=' + str(round(syn,4)))
            # 显示图形
            plt.show()
        
    return un_sys, un_en, syn, series_ls

def hamming_distance(seq1, seq2):
    # 确保两个序列长度相同
    if len(seq1) != len(seq2):
        raise ValueError("Both sequences must have the same length")
    
    # 计算汉明距离
    distance = 0
    for bit1, bit2 in zip(seq1, seq2):
        # 如果两个位不同，则距离加1
        if bit1 != bit2:
            distance += 1
    
    return distance

def dis_mean(series):
    dis = 0
    lens = len(series)
    nums = lens * (lens - 1) / 2
    for i in series:
        for j in series:
            dis += hamming_distance(i,j)
    dis /= nums
    return dis




# def permute_matrix_rows(original_order, new_order):
#     # 原始和新顺序的长度（应该是3）
#     n = len(original_order)
#     # 创建一个映射字典
#     mapping = []
    
#     # 生成所有可能的状态（000到111）
#     for state in product([0, 1], repeat=n):
#         # 将元组转换为二进制字符串
#         original_state = ''.join(map(str, state))
#         #new_state = ''.join(original_state[i] for i in [original_order.index(new_order[j]) for j in range(n)])
#         # 根据新顺序重新排列状态位
#         #indices = np.where(original_order == new_order[j] for j in range(n))[0][0]  # 注意：这里假设每个元素都是唯一的

#         # 然后你可以使用这些索引来重新排列状态位
#         new_state = ''.join(original_state[i] for i in [np.where(original_order == new_order[j])[0][0] for j in range(n)])
        
#         # 将原始状态映射到新状态
#         mapping.append(new_state)
    
#     return mapping

# def add_missing_elements(array1, array2):
#     # 将数组转换为集合
#     set1 = set(array1)
#     set2 = set(array2)
#     array22 = array2.copy()
#     # 找出在set1中但不在set2中的元素
#     missing_elements = set1 - set2
    
#     # 将缺失的元素添加到array2中
#     array22 = np.concatenate((np.array(list(missing_elements)), array22))
    
#     return array22

# def decimal_to_binary(decimal, min_length=1):
#     if min_length == 0:
#         return ''
#     if decimal == 0:
#         return "0" if min_length == 1 else "0".zfill(min_length)
#     binary = ""
#     while decimal > 0:
#         binary = str(decimal % 2) + binary
#         decimal = decimal // 2
#     # 使用 zfill 确保二进制字符串至少有 min_length 长度
#     return binary.zfill(min_length)

# def nei_comb(candidate_system, F, I, noise):
#     neigbors = candidate_system
#     tpm_list = []
#     all_in_list = []
#     for i in candidate_system:
#         neigbors = np.concatenate((neigbors, I[i]))
#         tpm, size, ins = tpm_one(F[i], I[i], i, noise=noise)
#         tpm_list.append(tpm)
#         all_in_list.append(ins)
#     #neigbors = np.unique(neigbors)
#     seen = set()
#     neigbors_un = [x for x in neigbors if not (x in seen or seen.add(x))]
# #     print(neigbors_un)
#     neigbors_ = np.setdiff1d(neigbors, candidate_system)
#     return neigbors_un, tpm_list, neigbors_, all_in_list

# def tpm_comb(candidate_system, F, I, noise):
#     neigbors, tpm_list, _, _ = nei_comb(candidate_system, F, I, noise)
#     matrix = np.ones([2**len(neigbors), 2**len(candidate_system)])
#     for i, s in enumerate(candidate_system):
#         new_arr = add_missing_elements(neigbors, I[s])
# #         print("new_arr")
# #         print(new_arr)
#         while tpm_list[i].shape[0] < matrix.shape[0]:
#             tpm_list[i] = np.tile(tpm_list[i], (2, 1))
#         mapping = permute_matrix_rows(neigbors, new_arr)
#         new_id = [int(index,2) for index in mapping]
#         tpm_list[i] = tpm_list[i][new_id]
#     for i in range(matrix.shape[1]):
#         ind_str = decimal_to_binary(i, min_length=len(candidate_system))
#         for k, j in enumerate(ind_str):
#             matrix[:, i] *= tpm_list[k][:, int(j)] 
#     return neigbors, matrix


# def tpm_ei(tpm, log_base = 2):
#     # marginal distribution of y given x ~ Unifrom Dist
#     puy = tpm.sum(axis=0)
#     n = tpm.shape[0]
#     # replace 0 to a small positive number to avoid log error
#     eps = 1E-10
#     tpm_e = np.where(tpm==0, eps, tpm)
#     puy_e = np.where(tpm==0, eps, puy)
    
#     # calculate EI of specific x
#     ei_x = (np.log2(n * tpm_e / puy_e) / np.log2(log_base)  * tpm).sum(axis=1)
    
#     # calculate total EI
#     ei_all = ei_x.mean()
#     return ei_all


# def synergy(markov_matrix, mech_size, en_size):
#     ei_all = condi_ei(markov_matrix, mech_size, en_size)
#     un = unique(markov_matrix, mech_size, en_size)[0]
#     syn = ei_all - un
#     return syn

# def condi_ei(markov_matrix, mech_size, en_size):
#     ei = 0
#     state_size = 2**mech_size
#     state_en = 2**en_size
#     for e in range(state_en):
#         local_markov = np.zeros([state_size, state_size])
#         for num in range(state_size):
#             binary_string = decimal_to_binary(num, min_length=mech_size)
#             padded_binary_string = binary_string + decimal_to_binary(e, min_length=en_size)
#             binary_array = [int(bit) for bit in padded_binary_string] 
#             pattern = int(''.join(str(cell) for cell in binary_array), 2)
#             local_markov[num, :] = markov_matrix[pattern, :]
#         ei += tpm_ei(local_markov)
#     ei = ei / state_en
#     return ei 

# def unique(markov_matrix, mech_size, en_size):
#     state_size = 2**mech_size
#     state_en = 2**en_size
#     mixed_markov = np.zeros([state_size, state_size])
#     for e in range(state_en):
#         local_markov = np.zeros([state_size, state_size])
#         for num in range(state_size):
#             binary_string = decimal_to_binary(num, min_length=mech_size)
#             padded_binary_string = binary_string + decimal_to_binary(e, min_length=en_size)
#             binary_array = [int(bit) for bit in padded_binary_string] 
#             pattern = int(''.join(str(cell) for cell in binary_array), 2)
#             local_markov[num, :] = markov_matrix[pattern, :]
#         mixed_markov += local_markov
#     ei = tpm_ei(mixed_markov / state_en)
#     return ei, mixed_markov

# def en_unique(markov_matrix, mech_size, en_size):
#     state_size = 2**mech_size
#     state_en = 2**en_size
#     mixed_markov = np.zeros([state_en, state_size])
#     for s in range(state_size):
#         local_markov = np.zeros([state_en, state_size])
#         for num in range(state_en):
#             binary_string = decimal_to_binary(num, min_length=en_size)
#             padded_binary_string = decimal_to_binary(s, min_length=state_size) + binary_string
#             binary_array = [int(bit) for bit in padded_binary_string] 
#             pattern = int(''.join(str(cell) for cell in binary_array), 2)
#             local_markov[num, :] = markov_matrix[pattern, :]
#         mixed_markov += local_markov
#     ei,det,deg,eff,det_c,deg_c = tpm_ei_new(mixed_markov / state_size)
#     return ei, mixed_markov

# def partition(tpm, inputs, eu_list):
#     in_rest = np.setdiff1d(inputs, eu_list)
# #     print(tpm)
# #     print(in_rest)
#     #res_ind = np.where(inputs==in_rest)[0]
#     #print(res_ind)
#     rows = 2**len(in_rest)
#     tpm_p = np.zeros([rows, 2])
#     for i in range(tpm.shape[0]):
#         strs = decimal_to_binary(i, min_length=len(inputs))
#         #print(np.where(inputs==in_rest[0]))
#         str_res = [strs[np.where(inputs==k)[0][0]] for k in in_rest]
#         str_res = ''.join(str_res)
# #         print(str_res)
#         num = int(str_res, 2)
#         tpm_p[num,:] += tpm[i, :]
# #         print(tpm_p)
#     tpm_p /= 2**len(eu_list)
#     return tpm_p

# def un_comb(candidate_system, F, I, noise):
#     neigbors, tpm_list, neighbors_, all_in_list = nei_comb(candidate_system, F, I, noise)
#     un = 0
#     set1 = set(neighbors_)
#     for i, s in enumerate(candidate_system):
# #         print(un)
#         set2 = set(I[s])
#         intersection = set1.intersection(set2)
#         eu_list = list(intersection)
#         if eu_list == []:
#             un += tpm_ei_new(tpm_list[i])[0]
#         else:
#             un += tpm_ei_new(partition(tpm_list[i], all_in_list[i], eu_list))[0]
#     return un

# def syn_comb(candidate_system, F, I, noise):
#     neigbors, tpm_list, neighbors_, all_in_list = nei_comb(candidate_system, F, I, noise)
#     syn = 0
#     set1 = set(neighbors_)
#     for i, s in enumerate(candidate_system):
#         inputs = all_in_list[i]
#         set2 = set(I[s])
#         intersection = set1.intersection(set2)
#         eu_list = list(intersection)
#         if len(eu_list) > 0:
#             condi = 0
#             in_rest = np.setdiff1d(inputs, eu_list)
#             rows = 2**len(in_rest)
#             for j in range(2**len(eu_list)):
#                 tpm_p = np.zeros([rows, 2])
#                 str_en = decimal_to_binary(j, min_length=len(eu_list))
#                 for k in range(2**len(inputs)):
#                     strs_all = decimal_to_binary(k, min_length=len(inputs))
#                     str_eu = [strs_all[np.where(inputs==e)[0][0]] for e in eu_list]
#                     str_eu = ''.join(str_eu)
#                     if str_eu == str_en:
#                         num = int(strs_all, 2)
#                         str_in = [strs_all[np.where(inputs==n)[0][0]] for n in in_rest]
#                         str_in = ''.join(str_in)
#                         num_in = int(str_in, 2)
# #                         print("sys:"+str(s))
# #                         print("tpm_p:"+str(tpm_p.shape))
# #                         print("num_in:"+str(num_in))
# #                         print("tpm_all:"+str(tpm_list[i].shape))
# #                         print("num:"+str(num))
#                         tpm_p[num_in, :] = tpm_list[i][num, :]
#                 condi += tpm_ei_new(tpm_p)[0]
#             syn += condi / 2**len(eu_list)
#             syn -= tpm_ei_new(partition(tpm_list[i], all_in_list[i], eu_list))[0]
#     return syn