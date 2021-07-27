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
P0 = [V0]

pi_v = 0.1955
k_v1 = 0.1955

model_args = (pi_v, k_v1)

#######   Viremia  log10 copias/ml ##########
dadosViremiaLog10 = pd.read_csv('../data/Viral_load.csv',',')

dias_de_simulação = 35

t=range(dias_de_simulação)    
y,d=integrate.odeint(immune_response_v1, P0, t, args=(model_args), full_output=1, printmessg=True)

plt.figure('Viremia')
plt.plot(t, dadosViremiaLog10['Viral_load'], 'o', label='data', linewidth=4)
plt.xlim(0.0,dias_de_simulação)
plt.plot(t,y[:,0],label='Curva gerada pelo modelo',linewidth=1.5, linestyle="-")
plt.xlabel('Tempo (dias)')
plt.ylabel('Viremia')
plt.legend()
#plt.savefig(output_path+'Viremia.pdf',bbox_inches='tight',dpi = 300)
plt.show()
