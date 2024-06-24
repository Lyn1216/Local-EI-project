import numpy as np
import func.entropy_estimators as ee

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