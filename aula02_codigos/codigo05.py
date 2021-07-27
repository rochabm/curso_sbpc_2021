#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from NovoModelo import *
from scipy.stats.mstats import gmean
import pylab
import math
sns.set_style("whitegrid")

output_path= './output/'


V0 = 1.5e0
Ap0 = 0.1
Apm0=0.0
Tkn0=0.0
Tke0=0.0
I0=0.0
C0=0.0
P0 = [V0,Ap0,Apm0,Tkn0,Tke0,I0,C0]

pi_v = 0.1955
k_v3 = 0.1
alpha_Ap = 1.87E-03
Ap_0 = 1.0e6
beta_Ap = 2.00E-03
k_ap1 = 0.8  
k_ap2 = 40.0
delta_Apm = 1.0e-04
k_v2= 9.5e-1
beta_Tk=2.10152618e-01
pi_T=1.431849023090428446e-05
k_te1=1.0E-08 
delta_te= 0.0003 
Tkn_0 = 5.0e5
k_apm= 5.36139617e-01
k_tk=2.10152618e-01
pi_c_apm=7.43773673e-00
pi_c_i= 1.97895565e-01
pi_c_tke=4.730172e-02
delta_c=8.26307952e-05

model_args = (pi_v, k_v3, alpha_Ap, Ap_0,beta_Ap,k_ap1,k_ap2,delta_Apm,k_v2, beta_Tk,pi_T,k_te1, delta_te,Tkn_0,k_apm, k_tk, pi_c_apm, pi_c_i, pi_c_tke, delta_c)

#######   Viremia  log10 copias/ml ##########
dadosViremiaLog10 = pd.read_csv('../data/Viral_load.csv',',')
dadosCitocinaObitos = pd.read_csv('../data/IL6_non-survivors_19.csv',',')
dadosCitocinaSobreviventes = pd.read_csv('../data/IL6_survivors_19.csv',',')


dias_de_simulação = 35

t=range(dias_de_simulação)    
y,d=integrate.odeint(immune_response_v5, P0, t, args=(model_args), full_output=1, printmessg=True)

plt.figure('Viremia')
plt.plot(t, dadosViremiaLog10['Viral_load'], 'o', label='data', linewidth=4)
plt.xlim(0.0,dias_de_simulação)
plt.plot(t,y[:,0],label='Curva gerada pelo modelo',linewidth=1.5, linestyle="-")
plt.xlabel('Tempo (dias)')
plt.ylabel('Viremia')
plt.legend()
#plt.savefig(output_path+'Viremia.pdf',bbox_inches='tight',dpi = 300)


plt.figure('Apcs')
plt.xlim(0.0,dias_de_simulação)
plt.plot(t,y[:,1],label='Curva gerada pelo modelo',linewidth=1.5, linestyle="-")
plt.xlabel('Tempo (dias)')
plt.ylabel('Apcs')
plt.legend()
#plt.savefig(output_path+'Viremia.pdf',bbox_inches='tight',dpi = 300)


plt.figure('C')
plt.plot(dadosCitocinaSobreviventes['Day']+5, dadosCitocinaSobreviventes['IL6(pg/mL)'], 'o', label='data', linewidth=4)
plt.xlim(0.0,dias_de_simulação)
#plt.ylim(0.0,8.0)
plt.plot(t,y[:,6],label='C',linewidth=1.5, linestyle="-")
plt.xlabel('Tempo pós-vacinação (dias)')
plt.ylabel('C')
plt.legend()
#plt.savefig(output_path+'C.pdf',bbox_inches='tight',dpi = 300)

plt.figure('C')
dadosCitocinaObitos.plot.scatter(x='Day',y='IL6(pg/mL)',color='y',label='Dados experimentais(Obito)')
plt.xlim(0.0,dias_de_simulação)
#plt.ylim(0.0,8.0)
plt.plot(t,y[:,6],label='C',linewidth=1.5, linestyle="-")
plt.xlabel('Tempo pós-vacinação (dias)')
plt.ylabel('C')
plt.legend()
#plt.savefig(output_path+'C2.pdf',bbox_inches='tight',dpi = 300)

plt.show()
