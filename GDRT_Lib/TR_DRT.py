import numpy as np
from scipy.optimize import minimize, Bounds
from .calculate_EIS import calc_A, omega

# TODO make class TR_DRT Where C and L (make constraints for L instead of fixing it) Compare results with already shown data and tikhonov minimization are instances
def TR_DRT(tau, Z_exp, x0=None, el=1e-2, method="SLSQP", lguess=1.3e-7, Cguess=0.15):
    gammaRC, R_inf, C = Tikhonov_minimization(tau, Z_exp, x0, el, method, lguess, Cguess)
    return gammaRC, R_inf, C

def Tikhonov_minimization(tau, Z_exp, x0, el, method, l, Cguess):
    N = len(Z_exp)
    # boundary conditions capacitive element
    bc = [1e-12, np.inf]
    Z_exp_re=np.real(Z_exp)
    Z_exp_im=np.imag(Z_exp)
    A_RC_re = calc_A(tau,option='RC_Re')
    A_RC_im = calc_A(tau,option='RC_Im')
    w = omega(tau)
    if x0 is None:
        x0 = initial_guess(tau, Z_exp, Cguess)
    bounds = Bounds(np.append(np.zeros(x0[:-1].shape),[bc[0]]), 
                    np.append(np.abs(Z_exp).max()*np.ones(x0[:-1].shape),[bc[1]]))
    result = minimize(SC, x0, args=(w, Z_exp_re, Z_exp_im, A_RC_re, A_RC_im, el, l), method=method,
                      bounds = bounds, options={'disp': False, 'ftol':1e-9,'maxiter':800})
    
    # print('Min result:', result.jac)
    gamma_R_inf = result.x
    R_inf = gamma_R_inf[N]
    # L = gamma_R_inf[N+2]
    C = gamma_R_inf[N+1]
    # n = gamma_R_inf[N+2]
    gammaRC = gamma_R_inf[:N]
    return gammaRC, R_inf, C

# Fixed inductance
def SC(gamma_R_inf, w, Z_exp_re, Z_exp_im, A_RC_re, A_RC_im, el, l):
    N = len(Z_exp_re)
    MSE_re = np.sum((gamma_R_inf[N] + np.matmul(A_RC_re, gamma_R_inf[:N]) - Z_exp_re)**2)
    MSE_im = np.sum(((-1/(gamma_R_inf[(N+1)]*w)) + l*w+ np.matmul(A_RC_im, gamma_R_inf[:N]) - Z_exp_im)**2)
    reg_term = el/2*np.sum(gamma_R_inf[:N]**2)
    obj = MSE_re + MSE_im + reg_term
    return obj

# # Fixed inductance and CPE
# def SCCPE(gamma_R_inf, w, Z_exp_re, Z_exp_im, A_RC_re, A_RC_im, el):
#     N = len(Z_exp_re)
#     MSE_re = np.sum((gamma_R_inf[N] + np.matmul(A_RC_re, gamma_R_inf[:N]) - Z_exp_re)**2)
#     MSE_im = np.sum(((-1/np.power(gamma_R_inf[(N+1)]*w,gamma_R_inf[(N+2)])) + 1.3e-7*w + np.matmul(A_RC_im, gamma_R_inf[:N]) - Z_exp_im)**2)
#     reg_term = el/2*np.sum(gamma_R_inf[:N]**2)
#     obj = MSE_re + MSE_im + reg_term
#     return obj

# def S(gamma_R_inf, w, Z_exp_re, Z_exp_im, A_RC_re, A_RC_im, el):
#     N = len(Z_exp_re)
#     MSE_re = np.sum((gamma_R_inf[N] + np.matmul(A_RC_re, gamma_R_inf[:N]) - Z_exp_re)**2)
#     MSE_im = np.sum((np.matmul(A_RC_im, gamma_R_inf[:N]) - Z_exp_im)**2)
#     reg_term = el/2*np.sum(gamma_R_inf[:N]**2)
#     obj = MSE_re + MSE_im + reg_term
#     return obj


# def SCL(gamma_R_inf, w, Z_exp_re, Z_exp_im, A_RC_re, A_RC_im, el):
#     N = len(Z_exp_re)
#     MSE_re = np.sum((gamma_R_inf[N] + np.matmul(A_RC_re, gamma_R_inf[:N]) - Z_exp_re)**2)
#     MSE_im = np.sum(((-1/(gamma_R_inf[(N+1)]*w)) + gamma_R_inf[(N+2)]*w + np.matmul(A_RC_im, gamma_R_inf[:N]) - Z_exp_im)**2)
#     reg_term = el/2*np.sum(gamma_R_inf[:N]**2)
#     obj = MSE_re + MSE_im + reg_term
#     return obj

def initial_guess(tau, Z_exp,Cguess):
    N=len(Z_exp)
    x0 = np.zeros(N+2)
    x0[N] = np.abs(Z_exp)[-1]
    x0[N+1] = Cguess
    # x0[N+2] = 0.88
    # x0[N+2] = 4e-7
    return x0

def errFit(hess_inv, resVariance):
    return np.sqrt( np.diag( hess_inv * resVariance))