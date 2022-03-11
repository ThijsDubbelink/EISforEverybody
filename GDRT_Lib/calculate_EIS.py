import numpy as np
from math import pi, log

''' Check paper Zhang et al. to improve this calculation'''
def calc_A(tau, option=False):
    omega = 2*pi/tau
    N_freqs = tau.size
    out_A_re = np.zeros((N_freqs, N_freqs))
    if option is False:
        print('choose an option')
    for p in range(0, N_freqs):
        for q in range(0, N_freqs):
            if option == 'RC_Re':
                Abas = -1/(1+(omega[p]*tau[q])**2)
            if option == 'RC_Im':
                Abas = (omega[p]*tau[q])/(1+(omega[p]*tau[q])**2)
            if q == 0:
                out_A_re[p, q] = 0.5*Abas*log(tau[q+1]/tau[q])
            elif q == N_freqs-1:
                out_A_re[p, q] = 0.5*Abas*log(tau[q]/tau[q-1])
            else:
                out_A_re[p, q] = 0.5*Abas*log(tau[q+1]/tau[q-1])
    return out_A_re

def omega(tau):
    return 2*pi/tau

def calculate_EIS(tau, gammaRC, R_inf, C, l):
    A_RC_re = calc_A(tau, option='RC_Re')
    A_RC_im = calc_A(tau, option='RC_Im')
    w = omega(tau)
    Z_cal = R_inf + np.matmul(A_RC_re, gammaRC) + 1j*(l*w-1/(w*C)+ np.matmul(A_RC_im, gammaRC))
    # Z_cal = R_inf + np.matmul(A_RC_re, gammaRC) + 1j*(np.matmul(A_RC_im, gammaRC))
    return Z_cal
