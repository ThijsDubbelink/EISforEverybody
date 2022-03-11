# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 21:29:48 2021
@author: Thijs Dubbelink
"""

import pandas as pd
import glob

'''ToDo change output to tau vector and Z vector'''
class Import:
             
    def __init__(self,path,ext='.csv'):
        self.path = path
        self.ext = ext
        self.flnms, self.LOF, self.DFL = [],[],[]
        self.lof()
        self.n = len(self.flnms)
        
    # def __init__(self,path,ext='.txt'):
    #     self.path = path
    #     self.ext = ext
    #     self.flnms, self.LOF = self.lof()
    #     self.n = len(self.flnms)
    #     self.EIS = []
    #     self.AA = 60.2
    
    def lof(self):
        flnmst = glob.glob(self.path+"/*"+self.ext)
        for file in flnmst:
            idb = file.find(self.path[-10:])+11
            ide = file.find(self.ext)
            self.LOF.append(file[idb:ide])
            self.flnms.append(file[:ide])
            
    def BiolAuto(self,freqlim=10e6):
        for fl in self.flnms:
            df = pd.read_csv(fl+self.ext, skiprows=1, delim_whitespace=True,
                             usecols=[0,1,2], names = ['freq','Z_re','Z_im'])
            df['-Z_im'] = df['Z_im']*-1
            del df['Z_im']
            # df[['Z_re','-Z_im']] *= self.AA
            df = df.iloc[::-1]
            df = df.loc[df['freq']<freqlim]
            self.DFL.append(df)
        
    def BiolCycAlt(self,freqlim=10e6):
        for fl in self.flnms:            
            df = pd.read_csv(fl+self.ext, delim_whitespace=True, skiprows=0)
            df = df.rename(columns={"Re(Z)/Ohm": "Z_re",
                                    "-Im(Z)/Ohm": "Z_im",
                                    "freq/Hz": "freq",
                                    "cycle number": "Cyc#"})
            df['-Z_im'] = df['Z_im']*-1
            # del df['Z_im']
            # df[['Z_re','-Z_im']] *= self.AA
            df = df.loc[(df['freq']<freqlim)&(df['freq']>1e-3)]
            self.DFL.append(df)

    def BiolCyc(self):
        for fl in self.flnms:
            # print(fl+self.ext)
            fl += self.ext
            with open(fl, 'r') as f:                                            #open file as f
                line = f.readlines()[1]
                idb = line.find(':')
                skiplines = int(line[idb+1:])
            df = pd.read_csv(fl, sep='\t' , skiprows=skiplines-1,
                             encoding="ISO-8859-1")
            self.DFL.append(df)
        
    def Filtercyceis(self,freqlim=10e6):
        self.BiolCyc()
        for df in self.DFL:
            df['type'] = ""
            df['EIS mode'] = ""
            df['EIS mode'].loc[(df['freq/Hz']<freqlim)&(df['freq/Hz']>1e-3)] = 'EIS'
            cur = df.loc[df['I/mA']>1e-3].groupby(['time/s'])['I/mA'].mean().abs().round(1).mode().values[0]            
            df['type'].loc[(df['I/mA']>0.95*cur)] = 'Charge'
            df['type'].loc[(df['I/mA']<-0.95*cur)] = 'Discharge'
            cyc = df['cycle number'].loc[df['type']=='Discharge'].drop_duplicates().values
            cycn = list(range(1,len(cyc)+1))
            df['cyc num new'] = float("nan")
            for c in range(len(cyc)):
                df['cyc num new'].loc[df['cycle number']==cyc[c]] = cycn[c]
            df['cyc num new'] = df['cyc num new'].fillna(method='bfill')
            # df = df.loc[df['cyc num new']!=float("nan")]
            df = df[df['cyc num new'].notna()]
''' drop na for all non finished EIS cycles'''
                