import numpy as np
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
    
def tpm_ei(tpm, log_base = 2):
    # marginal distribution of y given x ~ Unifrom Dist
    puy = tpm.sum(axis=0)
    n = tpm.shape[0]
    # replace 0 to a small positive number to avoid log error
    eps = 1E-10
    tpm_e = np.where(tpm==0, eps, tpm)
    puy_e = np.where(puy==0, eps, puy)
    
    # calculate EI of specific x
    ei_x = (np.log2(n * tpm_e / puy_e) / np.log2(log_base)  * tpm).sum(axis=1)
    
    # calculate total EI
    ei_all = ei_x.mean()
    return ei_all

# def tpm_ei_new(tpm, log_base = 2):
#     # marginal distribution of y given x ~ Unifrom Dist
#     puy = tpm.sum(axis=0)
#     n = tpm.shape[0]
#     # replace 0 to a small positive number to avoid log error
#     eps = 1E-10
#     tpm_e = np.where(tpm==0, eps, tpm)
#     puy_e = np.where(tpm==0, eps, puy)
    
#     # calculate EI of specific x
#     ei_x = (np.log2(n * tpm_e / puy_e) / np.log2(log_base)  * tpm).sum(axis=1)
    
#     det = np.log2(n) + (tpm * np.log2(tpm_e)).mean(axis=0)
#     deg = np.log2(n) + tpm.mean(axis=0) * np.log2(tpm_e.mean(axis=0))
#     # calculate total EI
#     ei_all = ei_x.mean()
#     return ei_all,det/np.log2(log_base),deg/np.log2(log_base)

def tpm_ei_new(tpm, log_base = 2):
    '''
    tpm: 输入的概率转移矩阵，可以是非方阵
    log_base：对数的底

    ei_all：EI的值
    det：EI中确定性的部分
    deg：EI中简并性的部分
    eff：有效性
    det_c：确定性系数
    deg_c：简并性系数
    '''
    # marginal distribution of y given x ~ Unifrom Dist
    puy = tpm.mean(axis=0)
    m,n = tpm.shape
    if m > n:
        q = n
    else:
        q = m
    
    # replace 0 to a positive number to avoid log error
    eps = 0.5
    tpm_e = np.where(tpm==0, eps, tpm)
    puy_e = np.where(puy==0, eps, puy)
    
    # calculate EI of specific x
    ei_x = (np.log2(tpm_e / puy_e) / np.log2(log_base)  * tpm).sum(axis=1)
    
    # calculate det and deg
    det = np.log2(n) + (tpm * np.log2(tpm_e)).sum(axis=1).mean(axis=0)
    deg = np.log2(n) + (puy * np.log2(puy_e)).sum()
    
    det = det / np.log2(log_base)
    deg = deg / np.log2(log_base)
    ei_all = ei_x.mean()
    if q > 1:
        det_c = det / np.log2(q) * np.log2(log_base)
        deg_c = deg / np.log2(q) * np.log2(log_base)
        eff = ei_all / np.log2(q) * np.log2(log_base)
    else:
        det_c = 0
        deg_c = 0
        eff = 0
    return ei_all,det,deg,eff,det_c,deg_c

def tpm_ei_new2(tpm, log_base = 2):
    # marginal distribution of y given x ~ Unifrom Dist
    puy = tpm.mean(axis=0)
    m,n = tpm.shape
    if m > n:
        q = n
    else:
        q = m
    
    # replace 0 to a positive number to avoid log error
    eps = 1e-5
    tpm_e = np.where(tpm==0, eps, tpm)
    puy_e = np.where(puy==0, eps, puy)
    
    # calculate EI of specific x
    ei_x = (np.log2(tpm_e / puy_e) / np.log2(log_base)  * tpm).sum(axis=1)
    
    # calculate det and deg
    det = (tpm * np.log2(tpm_e)).sum(axis=1).mean(axis=0)
    nondeg = (-puy * np.log2(puy_e)).sum()
    
    det = det / np.log2(log_base)
    nondeg = nondeg / np.log2(log_base)
    ei_all = ei_x.mean()
    
    return ei_all,det,nondeg

def condi_ei(markov_matrix,mech_size,state=1):
    ei = 0
    state_size = (state + 1)**(mech_size)
    for e1 in ['0','1']:
        for e2 in ['0','1']:
            local_markov = np.zeros([state_size, state_size])
            for num in range(state_size):
                binary_string = trans10_to_base(num, base = state +1, min_length = mech_size)
                padded_binary_string = e1 + binary_string + e2
                binary_array = [int(bit) for bit in padded_binary_string] 
                pattern = int(''.join(str(cell) for cell in binary_array), state+1)
                local_markov[num, :] = markov_matrix[pattern, :]
            ei0 = tpm_ei(local_markov, log_base = state+1)
            ei += ei0
    return ei / 4

def unique_ca(markov_matrix, mech_size, state=1):
    state_size = (state + 1)**(mech_size)
    mixed_markov = np.zeros([state_size, state_size])
    for e1 in ['0','1']:
        for e2 in ['0','1']:
            local_markov = np.zeros([state_size, state_size])
            for num in range(state_size):
                binary_string = trans10_to_base(num, base = state +1, min_length = mech_size)
                padded_binary_string = e1 + binary_string + e2
                binary_array = [int(bit) for bit in padded_binary_string] 
                pattern = int(''.join(str(cell) for cell in binary_array), state+1)
                local_markov[num, :] = markov_matrix[pattern, :]
            mixed_markov += local_markov
    ei = tpm_ei(mixed_markov / 4)
    return ei, mixed_markov

def en_unique_ca(markov_matrix, mech_size):
    state_size = 2**mech_size
    mixed_markov = np.zeros([4, state_size])
    for num in range(state_size):
        local_markov = np.zeros([4, state_size])
        for e1 in ['0','1']:
            for e2 in ['0','1']:
                binary_string = trans10_to_base(num, base = 2, min_length = mech_size)
                padded_binary_string = e1 + binary_string + e2 
                binary_array = [int(bit) for bit in padded_binary_string] 
                pattern = int(''.join(str(cell) for cell in binary_array), 2)
                local_markov[num, :] = markov_matrix[pattern, :]
        mixed_markov += local_markov
    ei,det,deg,eff,det_c,deg_c = tpm_ei_new(mixed_markov / state_size)
    return ei, mixed_markov

def synergy_ca(markov_matrix, mech_size):
    ei_all = condi_ei(markov_matrix,mech_size)
    un = unique_ca(markov_matrix,mech_size)[0]
    syn = ei_all - un
    return syn