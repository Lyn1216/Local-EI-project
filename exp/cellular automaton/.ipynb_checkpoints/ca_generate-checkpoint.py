import sys
sys.path.append("../..")
import numpy as np
import random
import matplotlib.pyplot as plt
import func.entropy_estimators as ee

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

def generate_markov(p0, rule, middle_size = 1):
    markov_matrix = np.zeros([2**(middle_size+2),2**middle_size])
    binary_string = bin(rule)[2:]  # 将十进制数转换为二进制字符串，并去掉前缀'0b'
    padding_length = 8 - len(binary_string)
    padded_binary_string = '0' * padding_length + binary_string
    
    markov_matrix1 = generate_markov_middleone(p0,rule)
    #print(markov_matrix1)
    
    for i in range(2**(middle_size+2)):
        ii = trans10_to_base(i,min_length=middle_size+2)
        for j in range(2**middle_size):
            j_str = trans10_to_base(j,base=2,min_length=middle_size)
            j_num = [int(q) for q in j_str]
            element = 1
            for k,kk in enumerate(j_num):
                jj = ii[k:k+3]
                jj_indice = int(''.join(str(cell) for cell in jj), 2)
                element = element * markov_matrix1[jj_indice,kk]
       
            markov_matrix[i,j] = element
            #markov_matrix[i,:] += p0 / 2**middle_size
    return markov_matrix