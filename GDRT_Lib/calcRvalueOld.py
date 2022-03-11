# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 21:37:21 2020
'''Why doesnt this function not calculate peak at the lowest frequency?'''
@author: Thijs Dubbelink
"""
import numpy as np
import pandas as pd

def calcRvalue(arr, tau, show=False):
    peakidx = findLocalMin(arr)
    peakmx = findLocalMax(arr)
    rdrt = [[],[],[],[],[]]
    for x in range(1, len(peakidx)):
        gr = arr[peakidx[x-1]:peakidx[x]]
        fr = np.diff(np.log(tau[peakidx[x-1]:peakidx[x]]))
        taumx = tau[peakmx[x-1]]
        rdrt[0].append(np.sum(gr[1:]*-fr))
        rdrt[1].append(np.log10(tau[peakidx[x]]))
        rdrt[2].append(np.log10(tau[peakidx[x-1]]))
        rdrt[3].append(np.log10(taumx))
        rdrt[4].append(arr[peakmx[x-1]])
    rdf = pd.DataFrame({'rval': rdrt[0], 'tauminlow': rdrt[1], 'tauminhigh': rdrt[2], 'taumax': rdrt[3],'gammapeak': rdrt[4]})
    if show is True:
        print(rdf)
    return rdf
    
def findLocalMin(arr):  
  
    # Empty lists to store points of  
    # local minima  
    mn = []  
    if(len(arr)>1):
        if(arr[0] < arr[1]):  
            mn.append(0)  
      
        # Iterating over all points to check  
        # local maxima and local minima  
        for i in range(1, len(arr)-1):  
      
            # Condition for local minima  
            if(arr[i-1] > arr[i] < arr[i + 1]):  
                mn.append(i)  
        if(arr[-1] < arr[-2]):  
            mn.append(len(arr)-1)  
        
    return mn

def findLocalMax(arr):  
  
    # Empty lists to store points of  
    # local maxima
    mx = []  
    if(arr[0] > arr[1]):  
        mx.append(len(arr)-1)  
        
    # Iterating over all points to check  
    # local maxima and local minima  
    for i in range(1, len(arr)-1):  
  
        # Condition for local maxima  
        if(arr[i-1] < arr[i] > arr[i + 1]):  
            mx.append(i)  
            
    if(arr[-1] > arr[-2]):  
        mx.append(len(arr)-1)  
  
    return mx


