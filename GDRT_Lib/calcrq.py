# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 12:51:44 2021

'''Why doesnt this function not calculate peak at the lowest frequency?'''
@author: Thijs Dubbelink
"""
import numpy as np
import pandas as pd
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm 
from scipy.optimize import curve_fit

def calcrq(taux, arr, pksp, gn= 0.7, show=False, pth = ''):
    logtau = np.log10(taux)
    taux = taux[logtau>pksp[1]]
    print(min(taux))
    arr = arr[logtau>pksp[1]]
    logtau = np.log10(taux)
    r, q, n, tau, guess, lb, ub= [], [], [], [], [], [], []
    for x in range(1, len(pksp)-1):
        gtau = 10**logtau[(logtau>pksp[x])&(logtau<pksp[x+1])].mean()
        gr = arr[(logtau>pksp[x])&(logtau<pksp[x+1])].max()
        guess[len(guess):] = [gr, gtau, gn]
        lb[len(lb):] = [0,10**pksp[x],0]
        ub[len(ub):] = [np.inf,10**pksp[x+1],1]
    try:
        popt, pcov = curve_fit(multi_rq, taux, arr, guess, bounds=(lb,ub))
        
    except RuntimeError:
        print('Runtime Error in calcrq')
        popt = [float('nan')]*3*len(guess)
        pcov = [[float('nan')]*3*len(guess)]*len(guess)
        
    perr = np.sqrt(np.diag(pcov))
    rdf = pd.DataFrame()
    for x in range(1, len(pksp)-1):
        cr = popt[(x-1)*3]
        crcov = perr[(x-1)*3]
        ctau = popt[(x-1)*3+1]
        ctaucov = perr[(x-1)*3+1]
        cn = popt[(x-1)*3+2]
        cncov = perr[(x-1)*3+2]
        
        '''Is this correct'''
        cq = ctau**cn/cr
        rdf['tau_'+str(x)] = [ctau]
        rdf['relerr_tau_'+str(x)] = [ctaucov/ctau*100]
        rdf['R_'+str(x)] = [cr]
        rdf['relerr_R_'+str(x)] = [crcov/cr*100]
        rdf['CPE_'+str(x)] = [cq]
        rdf['n_'+str(x)] = [cn]
        rdf['relerr_n_'+str(x)] = [cncov/cn*100]
        
    if show == 'Dataframe':
        print(rdf)
    if show == 'Plot':
        fig , ax = plt.subplots(1,1, figsize=(7,6))
        ax.semilogx(taux, arr, 'o', label='Data', color = 'black')
        fitx = np.power(10,np.linspace(pksp[1],pksp[-1],1000)) 
        for x in range(1, len(pksp)-1):
            ax.semilogx(fitx, rq(fitx, popt[(x-1)*3],popt[(x-1)*3+1],popt[(x-1)*3+2]), 
                     color = cm.hsv((x-1)/(len(pksp)-2)),
                     linestyle = 'solid', linewidth=2, 
                     label='R%s = %.4f $[\u03A9*cm^2]$' %((x-1),rdf['R_'+str(x)]))
        ax.semilogx(fitx, multi_rq(fitx, *popt), color = 'black', 
                linestyle = 'dashed', linewidth=2, 
                label = 'Tot R')
        ax.set_xlabel('$\u03C4$ $[s]$', fontsize = 15)
        ax.set_ylabel('$\u03B3$ $[\u03A9*cm^2]$', fontsize = 15)
        fig.legend()
        # fig.show()
        if pth != '':
            fig.savefig(pth+".png")
    return rdf
    
def rq(x, R, TAU, phi):
    return R*np.sin(phi*np.pi)/(np.pi*2)/(np.cosh(phi*np.log(TAU/x))+np.cos(phi*np.pi))

def multi_rq(x, *pars):
    gs = []
    for xj in range(int(len(pars)/3)):
        tx = xj*3
        g = rq(x, pars[tx], pars[tx+1], pars[tx+2])
        gs.append(g)
    return sum(gs)



