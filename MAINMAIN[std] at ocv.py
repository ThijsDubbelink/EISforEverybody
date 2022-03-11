# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 16:04:32 2021

@author: Thijs Dubbelink
"""
import os



import pandas as pd
import numpy as np
from math import pi, sin, cos
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import GDRT_Lib
import sqlite3
import matplotlib.cm as cm 
import math

dupl = 2
pth = r'' # filelocation of .mpt data
savepth = pth
obj = GDRT_Lib.Import(pth, ext='.mpt')
obj.Filtercyceis(2e5)

obj.AA = np.ones(obj.n)*30.1
# %%

peaksep = [-np.inf,-3.8,-2.5,-1.4,-0.5,0.5,2]
rldf = []
plts = []
fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(18,9))
dfexp = pd.DataFrame()
for bat in range(obj.n):
    lg = 1e-7 * obj.AA[bat]
    cg = 0.2

    dfbat = pd.DataFrame()
    dfres = pd.DataFrame()
    batdf = obj.DFL[bat]
    batdf = batdf[batdf['cyc num new'].notna()]
    batdfc = batdf.loc[(batdf['EIS mode']=='EIS')]
    
    strt = batdfc.loc[batdfc['freq/Hz'].max()==batdfc['freq/Hz']].index.values
    nd = batdfc.loc[batdfc['freq/Hz'].min()==batdfc['freq/Hz']].index.values
    if len(strt)>len(nd):
        strt = strt[:len(nd)]
    for x in range(len(strt)):
        while nd[x]<strt[x]:
            nd = np.delete(nd,x)
    dfbat['BatsubID'] = len(strt)*[bat]
    # dfbat['idx'] = list(range(0,len(strt)))
    dfbat['BatID'] = len(strt)*[obj.LOF[bat]]
    dfbat['Volts'] = batdfc['Ewe/V'].loc[strt].values
    dfbat['cap'] = batdfc['Capacity/mA.h'].loc[strt].values
    dfbat['cyc'] = batdfc['cyc num new'].loc[strt].values
    rpts = len(dfbat.loc[dfbat['cyc']==dfbat['cyc'].min()])
    dfbat['idx'] = np.tile(np.linspace(0,rpts-1,rpts),len(strt)//rpts)
    # for n in range(7,len(strt),rpts): 
    for n in range(len(strt)):
        new = batdfc.loc[strt[n]:nd[n]]
        temp = new.iloc[::-1]
        tau = 1/temp['freq/Hz'].values
        Z_exp = temp['Re(Z)/Ohm'].values-1j*temp['-Im(Z)/Ohm'].values
        Z_exp*=obj.AA[bat]
        # Calculate distribution DRT
        gammaRC, R_inf, C = GDRT_Lib.TR_DRT(tau, Z_exp,lguess=lg, Cguess=cg, el=5e-3, method='SLSQP')
        # print('Rinf: ', R_inf, 'C: ', C)
        
        # filter frequency range
        Z_cal = GDRT_Lib.calculate_EIS(tau, gammaRC, R_inf, C, lg)
        '''This needs to be a function'''
        ss_res = np.sum((abs(Z_exp) - abs(Z_cal)) ** 2)
        ss_tot = np.sum((abs(Z_exp) - abs(np.mean(Z_exp))) ** 2)
        r2 = 1 - (ss_res / ss_tot)
        ax1.plot(np.real(Z_exp), -np.imag(Z_exp), "o", 
                  color = cm.tab20(bat//dupl/obj.n*dupl))
        ax1.plot(np.real(Z_cal), -np.imag(Z_cal), "x", color = 'black')
        ax1.legend(frameon=False, fontsize = 12)
        ax1.set_xlabel('$ReZ [\u03A9*cm^2]$', fontsize = 15)
        ax1.set_ylabel('$-ImZ [\u03A9*cm^2]$', fontsize = 15)
        ax1.axis("equal")
        ax2.plot(np.log10(tau), gammaRC, '.', linestyle = 'solid', 
                  color = cm.tab20(bat//dupl/obj.n*dupl), label = obj.LOF[bat])
        for pks in peaksep[1:]:
            ax2.axvline(pks, color='black')
        ax2.set_xlabel('$^{10}$log($\u03C4$) $[-]$', fontsize = 15)
        ax2.set_ylabel('$\u03B3$ $[\u03A9*cm^2]$', fontsize = 15)
        fig.legend(loc="upper right", bbox_to_anchor=(1,1), 
                borderaxespad=0, shadow=True, title='ID')
        # plt.show()

        # calculate resistances and capacitances
        newpath = pth+'/'+obj.LOF[bat]
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        rdf = GDRT_Lib.calcrq(tau, gammaRC, peaksep,
                              show = 'Plot',pth = newpath+'/Cyc='+str(dfbat['cyc'][n])+',V='+str(dfbat['Volts'][n]))
        
        
        rdf['R_0'] = R_inf
        rdf['R_-1'] = 0
        for rv in range(len(peaksep)-1):
            rdf['R_-1'] += rdf['R_'+str(rv)] 
        rdf['Rsquared'] = r2
        # print('Rsquared:',r2)
        dfres = pd.concat([dfres,rdf],ignore_index=True)
    dfbat = pd.concat([dfbat,dfres],axis=1)
    
    fig.savefig(savepth+'/Qualitative.png')
        
    dfexp = pd.concat([dfexp,dfbat],ignore_index=True)
dfexp.to_csv(pth+'/EISAnalysis.csv')
        # rinfdf = pd.DataFrame([[0, R_inf, 0, 0]], columns=rdf.columns)
        # rdf = pd.concat([rinfdf,rdf])
        # rdf.to_csv(pth+'/'+obj.LOF[x]+'rqR.csv')
        # rldf.append(rdf)   
    
    
    # # %%
    
    # fig1, ax = plt.subplots(1, len(peaksep)-1, figsize=(30,6))
    # fig2, ax4 = plt.subplots(1, 1, figsize=(6,6))
    # mx = 0
    # cl = ['red','black']
    # for j in range(obj.n):
    #     for x in range(len(peaksep)-1):
    #         r = rldf[j]['R'].loc[(rldf[j]['tau']>=10**peaksep[x]) &
    #                                   (rldf[j]['tau']<=10**peaksep[x+1])].values 
    #         r = sum(r)
    #         if r == 0:
    #             r = np.nan
    #         if r > mx:
    #             mx = r
            
    #         ax[x].plot(ml[j],r,'o', color = cl[j%dupl])
    #         ax[x].set_xlabel('Cycles (#)')
    #         ax[x].set_title('10$^{%.2f}$ < $\u03C4_{%s}$  < 10$^{%.2f}$ [s]'%(peaksep[x],x,peaksep[x+1]))       
    #     ax4.plot(ml[j],rldf[j]['R'].sum(), 'o', color = cl[j%dupl])
    # ax[0].set_ylabel('R $[\u03A9*cm^2]$')
    # ax4.set_ylim(0)
    
    # for xt in range(len(peaksep)-1):
    #     ax[xt].set_ylim(0,mx*1.1)
    # ax4.set_ylabel('Rtot $[\u03A9*cm^2]$')
    # ax4.set_xlabel('Cycles (#)')
    # fig1.savefig(savepth+'/QauntR_V='+str(np.mean(volts))+'.png')
    # fig2.savefig(savepth+'/QauntRtot_V='+str(np.mean(volts))+'.png')
    
    # # %%
    # clmns = []
    # df = pd.DataFrame()
    # df['Filename']=obj.LOF
    # for x2 in range(len(peaksep)-1):        
    #     r2 = []   
    #     for j2 in range(obj.n):
    #         r = rldf[j2]['R'].loc[(rldf[j2]['tau']>=peaksep[x2]) &
    #                                   (rldf[j2]['tau']<=peaksep[x2+1])].values
            
    #         r = sum(r)
    #         if r == 0:
    #             r = np.nan
    #         r2.append(r)
    #     df[' %.2f < logtau < %.2f [s]'%(peaksep[x2],peaksep[x2+1])] = r2
    
    # df.to_csv(pth+'/rqR.csv')
