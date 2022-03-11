# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 15:01:51 2021

@author: Labtop10
"""
import numpy as np
import matplotlib.pyplot as plt

def rq(x, R, TAU, phi):
    # return R/(2*np.pi)*np.sin((1-phi)*np.pi)*(1/(np.cosh(phi*np.log(TAU/x))-np.cos(1-phi)*np.pi))
    return np.sin(phi*np.pi)/(np.pi*2)/(np.cosh(phi*np.log(TAU/x))+np.cos(phi*np.pi))

# gammatemp = R[nb]/(2.*pi)*sin((1.-n[nb])*pi)*(1/(np.cosh(n[nb]*np.log(tau_plot/tautemp))-cos((1.-n[nb])*pi)))
fitx = np.power(10,np.linspace(2,-6,1000)) 
fig, ax = plt.subplots(1,1, figsize=(7,6))

ax.semilogx(fitx, rq(fitx,3.36,0.0012,0.7))