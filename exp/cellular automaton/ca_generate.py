import sys
sys.path.append("../..")
import numpy as np
import random
import matplotlib.pyplot as plt
import func.entropy_estimators as ee
from func.EI_calculation import unique_ca, en_unique_ca, synergy_ca

def trans10_to_base(number, base = 2, min_length=0):
    
    if min_length < 1:
        raise ValueError("Minimum length must be at least 1")
    
    if number == 0:
        return '0' * min_length

    digits = []
    while number > 0:
        digits.insert(0, str(number % base))
        number //= base

    # 将数字列表转换为字符串，并在前面填充字符以达到最小长度
    padded_digits = ''.join(digits).zfill(min_length)

    return padded_digits

def binary_vector_to_decimal(binary_vector):
    decimal_number = 0
    for index, bit in enumerate(binary_vector):
        decimal_number += bit * (2 ** index)
    return int(decimal_number)

def generate_markov_middleone(p0,rule):
    markov_matrix = np.zeros([8,2])
    binary_string = bin(rule)[2:]  # 将十进制数转换为二进制字符串，并去掉前缀'0b'
    padding_length = 8 - len(binary_string)
    padded_binary_string = '0' * padding_length + binary_string
    for i in range(8):
        ii = trans10_to_base(i,min_length=3)
        jj1 = int(''.join(str(cell) for cell in ii), 2)
        indice = int(padded_binary_string[jj1])
        markov_matrix[i,indice] = 1 - p0
        markov_matrix[i,1-indice] += p0
    return markov_matrix

def generate_markov(p0, rule, mech_size = 1):
    markov_matrix = np.zeros([2**(mech_size+2),2**mech_size])
    binary_string = bin(rule)[2:]  # 将十进制数转换为二进制字符串，并去掉前缀'0b'
    padding_length = 8 - len(binary_string)
    padded_binary_string = '0' * padding_length + binary_string
    
    markov_matrix1 = generate_markov_middleone(p0,rule)
    #print(markov_matrix1)
    
    for i in range(2**(mech_size+2)):
        ii = trans10_to_base(i,min_length=mech_size+2)
        for j in range(2**mech_size):
            j_str = trans10_to_base(j,base=2, min_length=mech_size)
            j_num = [int(q) for q in j_str]
            element = 1
            for k,kk in enumerate(j_num):
                jj = ii[k:k+3]
                jj_indice = int(''.join(str(cell) for cell in jj), 2)
                element = element * markov_matrix1[jj_indice,kk]
       
            markov_matrix[i,j] = element
            #markov_matrix[i,:] += p0 / 2**middle_size
    return markov_matrix

def move(strip, markov_matrix,state):
    grid1=np.zeros([1,len(strip)-2])
    cumulative_prob = np.cumsum(markov_matrix, axis=1)
    pattern = int(''.join(str(int(cell)) for cell in strip), state+1)
    value = random.random()
    indices = cumulative_prob[pattern,:] <= value
    interval_indices = np.argmin(indices)
    interval = trans10_to_base(interval_indices, base = state+1, min_length=len(strip)-2)
    grid1[0,:] = np.array(list(map(int, list(interval))))
    return grid1

def cellular_automaton_homo(rule, generations=120, size=120, p0_list=[0], mech_size=1, init=0, figure_show=False, vivid=False):
    period=len(p0_list)
    sub_size= size//period
#     # 初始化元胞状态
    if init == "random":
        current_generation = [random.randint(0, 1) for _ in range(size)]
    else:
        current_str = trans10_to_base(init, min_length=size)
        current_generation = [int(i) for i in current_str]
        
    showmatrix = np.zeros([generations+1,size])
    showmatrix[0,:] = current_generation
    
    markov_list = []
    for p0 in p0_list:
        markov_list.append(generate_markov(p0, rule, mech_size))
        
    un_mat = np.zeros([generations,size])
    syn_mat = np.zeros([generations,size])
    un_en_mat = np.zeros([generations,size])
    
    binary_string = bin(rule)[2:]  # 将十进制数转换为二进制字符串，并去掉前缀'0b'
    padding_length = 8 - len(binary_string)
    padded_binary_string = '0' * padding_length + binary_string
        

    for k in range(generations):
        next_generation = np.zeros(size)
        current_part = len(current_generation) // mech_size
        # 更新每个元胞的状态
        for i in range(size):
            index = i // sub_size
            markov_m = markov_list[index]
            if i > 0 and i < size-1 and vivid:
                un = unique_ca(markov_m, mech_size)[0]
                syn = synergy_ca(markov_m, mech_size)
                un_en = en_unique_ca(markov_m, mech_size)[0]
                un_mat[k,i] = un
                syn_mat[k,i] = syn
                un_en_mat[k,i] = un_en
            if i % mech_size == 0:
                if len(current_generation[i:i+mech_size+2])== mech_size+2:
                    next_generation[i+1:i+mech_size+1] = move(current_generation[i:i+mech_size+2],markov_m,state=1)
        current_generation = next_generation
        showmatrix[k+1,:] = current_generation
           
    
    if figure_show:
        plt.figure(figsize=(10,10)) 
        plt.imshow(showmatrix, aspect='auto')
        plt.show()
        plt.figure(figsize=(10,10)) 
        plt.imshow(un_mat[:, 1:size-mech_size], aspect='auto',cmap='hot')
        plt.colorbar()
        plt.figure(figsize=(10,10)) 
        plt.imshow(un_en_mat[:, 1:size-mech_size], aspect='auto',cmap='hot')
        plt.colorbar()
        plt.figure(figsize=(10,10)) 
        plt.imshow(syn_mat[:, 1:size-mech_size], aspect='auto',cmap='hot')
        plt.colorbar()
        plt.show()
    
    return showmatrix, un_mat, un_en_mat, syn_mat


# def cellular_automaton2(rule, generations = 100,size=100,p0_list=[0],middle_size = 1):
#     period=len(p0_list)
#     sub_size= size//period
#     # 初始化元胞状态
#     current_generation = [random.randint(0, 1) for _ in range(size)]
#     showmatrix = np.zeros([generations+1,size])
#     showmatrix[0,:] = current_generation
    
#     markov_list = []
#     for p0 in p0_list:
#         markov_list.append(generate_markov(p0,rule,middle_size))
        
#     ei_matrix=np.zeros([generations,size])
    
#     binary_string = bin(rule)[2:]  # 将十进制数转换为二进制字符串，并去掉前缀'0b'
#     padding_length = 8 - len(binary_string)
#     padded_binary_string = '0' * padding_length + binary_string
        

#     for k in range(generations):
#         next_generation = np.zeros(size)
#         current_part = len(current_generation) // middle_size
#         # 更新每个元胞的状态
#         for i in range(size):
#             index = i // sub_size
#             markov_m = markov_list[index]
#             if i > 0 and i < size - middle_size:
#                 e1 = str(int(current_generation[i - 1]))
#                 e2 = str(int(current_generation[i + middle_size]))
#                 ei, _ = local_ei_cell(middle_size, e1, e2, markov_m, state=1)
#                 ei_matrix[k,i] = ei 
#             if i % middle_size == 0:
#                 if len(current_generation[i: i+middle_size+2]) == middle_size+2:
#                     next_generation[i+1: i+middle_size+1] = move(current_generation[i:i+middle_size+2],markov_m,state=1)
#         current_generation = next_generation
#         showmatrix[k+1,:] = current_generation
    
#     plt.figure(figsize=(10,10)) 
#     plt.imshow(showmatrix, aspect='auto')
#     plt.show()
#     plt.figure(figsize=(10,10)) 
#     plt.imshow(ei_matrix[:, 1:size-middle_size], aspect='auto',cmap='hot')
#     plt.colorbar()
#     plt.show()
    
#     return showmatrix,ei_matrix


