# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 21:37:21 2020
'''Why doesnt this function not calculate peak at the lowest frequency?'''
@author: Thijs Dubbelink
"""
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm 
from scipy.optimize import curve_fit

Rcutoff = 1e-5
def calcRvalue(tau, arr, R_inf, show=False, pth = ''):
    logtau = np.log10(tau)
    peakmx = findLocalMax(arr)
    rdrt = [[-np.inf],[np.nan],[R_inf]]
    guess = []
    lb = []
    ub = []
    gsig = 0.4
    for x in range(0, len(peakmx)):
        mn = logtau[peakmx[x]]
        amp = arr[peakmx[x]]
        rdrt[0].append(mn)
        rdrt[1].append(amp)
        guess[len(guess):] = [amp, mn, gsig]
        lb[len(lb):] = [0,-np.inf,0]
        ub[len(ub):] = [np.inf,np.inf,np.inf]
    popt, pcov = curve_fit(multi_gaussian, logtau, arr, guess, bounds=(lb,ub))
    
    for x in range(0, len(peakmx)):
        r = areagaus(popt[x*3],popt[x*3+1],popt[x*3+2])
        rdrt[2].append(r)           
    rdf = pd.DataFrame({'taumax': rdrt[0],'gammapeak': rdrt[1],
                        'GaussR': rdrt[2]})
    rdf['Cap'] = 10**rdf['taumax'].div(rdf['GaussR'])
    if show == 'Dataframe':
        print(rdf)
    if show == 'Plot':
        fig , ax = plt.subplots(1,1, figsize=(7,6))
        ax.semilogx(tau, arr, 'o', label='Data, R_inf %.4f $[\u03A9]$' %rdf['GaussR'].loc[0], color = 'black')
        fitx = np.linspace(2,-6,1000) 
        for x in range(0, len(peakmx)):
            ax.semilogx(np.power(10,fitx), gaussian(fitx, popt[x*3],popt[x*3+1],popt[x*3+2]), 
                     color = cm.hsv(x/len(rdf)),
                     linestyle = 'solid', linewidth=2, 
                     label='R%s = %.4f $[\u03A9]$' %(x,rdf['GaussR'].loc[x+1]))
        ax.plot(np.power(10,fitx), multi_gaussian(fitx, *popt), color = 'black', 
                linestyle = 'dashed', linewidth=2, 
                label = 'Tot R = %.4f $[\u03A9]$' %rdf['GaussR'].sum())
        ax.set_xlabel('$\u03C4$ $[s]$', fontsize = 15)
        ax.set_ylabel('$\u03B3$ $[\u03A9]$', fontsize = 15)
        fig.legend()
        fig.show()
        if pth != '':
            fig.savefig(pth+".png")
    return rdf
    
def gaussian(x, A, x0, sig):
    return A*np.exp(-(x-x0)**2/(2*sig**2))

def areagaus(A, x0, sig):
    return A*math.sqrt(2*math.pi)*abs(sig)

def multi_gaussian(x, *pars):
    gs = []
    for xj in range(int(len(pars)/3)):
        tx = xj*3
        g = gaussian(x, pars[tx], pars[tx+1], pars[tx+2])
        gs.append(g)
    return sum(gs)

def findLocalMax(arr):  
  
    # Empty lists to store points of  
    # local maxima
    mx = []  
    if(arr[0] > arr[1]) & (arr[0]>Rcutoff):  
        mx.append(0)  
        
    # Iterating over all points to check  
    # local maxima and local minima  
    for i in range(1, len(arr)-1):  
  
        # Condition for local maxima  
        if(arr[i-1] < arr[i] > arr[i + 1]) & (arr[i]>Rcutoff):  
            mx.append(i)  
            
    if(arr[-1] > arr[-2]) & (arr[-1]>Rcutoff):  
        mx.append(len(arr)-1)  
  
    return mx

