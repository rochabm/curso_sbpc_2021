# Exemplo de ajuste de parametro de modelo SIR para COVID-19
# Codigo adaptado da seguinte referencia:
# https://medium.com/analytics-vidhya/coronavirus-in-italy-ode-model-an-parameter-optimization-forecast-with-python-c1769cf7a511
# Bernardo M. Rocha
# 13/05/2021

import sys
import numpy as np
import random
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.dates as mdates
from matplotlib import dates
from sklearn.metrics import mean_squared_error, r2_score
from datetime import datetime
from lmfit import minimize, Parameters, Parameter, report_fit

def deriv(y, t, N, ps):
    S, I, R = y
    try:
        beta_i = ps['beta_i'].value
        tau = ps['tau'].value
        gamma = ps['gamma'].value
    except:
        beta_i, beta_l, tau, gamma = ps
    
    beta = beta_i*(1.1-tau*t)
    #beta = beta_i
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

def odesol(y,t,N,ps):
    I0 = ps['i0'].value
    y0 = S0, I0, R0
    x = odeint(deriv, y0, t, args=(N, ps))
    return x
    
def residual(ps, ts, data):
    model = pd.DataFrame(odesol(y0,t,N,ps), columns=['S','I','R'])
    return (model['I'].values - data).ravel()

# dados reportados da italia
url = 'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv'
dfr = pd.read_csv(url,error_bad_lines=False)
dfr['data'] = pd.to_datetime(dfr['data'])
dfr['Days'] = (dfr.data.diff().dt.days+dfr.data.diff().dt.seconds/86400).cumsum().fillna(0).round(0)

# Total population, N.
N = 10e6

# Initial number of infected and recovered individuals, I0 and R0.
print(dfr.columns)
I0, R0 = dfr['totale_positivi'].min(), 0

# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0

ndays = 50

# Initial conditions vector
y0 = S0, I0, R0
#t = np.linspace(0, dfr.shape[0]-1, dfr.shape[0])
t = np.linspace(0, ndays, ndays) #dfr.shape[0])

# set parameters incluing bounds
params = Parameters()
params.add('i0', value=220, min=200, max=400)
params.add('beta_i', value= 0.37, min=0.1, max=0.5)
params.add('gamma', value= 0.11, min=0.1, max=0.2)
params.add('tau', value= 0.021, min=0.01, max=0.03)

#real data
data = dfr['totale_positivi'].values
data = data[:ndays]
print(np.shape(data))
print(data)

# fit model and find predicted values
result = minimize(residual, params, args=(t, data), method='leastsq')
final = data + result.residual.reshape(data.shape)

# display fitted statistics
report_fit(result)

# plot data and fitted curves
plt.plot(t, data, 'o',c='k', label='infectados (dados)')
plt.plot(t, final, '--', linewidth=2, c='red', label='infectados (modelo) best-fit-SIR');
#plt.yscale('log')
plt.xlabel('Dias')
plt.ylabel('Numero de Infectados')
plt.show()


sys.exit(0)

# ------------------------------------------------------------------------------------



### new
#fc_series = pd.Series(plotdata['Istd'])
#lower_series = pd.Series(plotdata['Ilow'])
#upper_series = pd.Series(plotdata['Isup'])

fig, ax = plt.subplots(1,2,figsize=(16,6))
fig.autofmt_xdate()

i=0
for axx in ax:
    sns.lineplot(x='data',y='totale_positivi', data=dfr, label='Total act. positives', lw=2, marker="o", size=50,legend='full',ax=ax[i])
    
    #ax[i].plot(fc_series, label='forecast')    
    #ax[i].fill_between(lower_series.index, lower_series, upper_series, color='k', alpha=.15)

    ax[i].set_ylabel('Units', fontsize=14)
    ax[i].set_xlabel('Date', fontsize=14)
    ax[i].xaxis.set_major_locator(mdates.AutoDateLocator())
    #ax[i].set(xticks=pd.date_range(plotdata.date.min(), periods=fc_series.shape[0]/7, freq='W'))
    ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%m.%d'))
    ax[i].set_ylim(0,100000)
    if i==1:
        ax[i].set_yscale('log')
        ax[i].set_ylim(100,200000)
        ax[i].set_ylabel('Log Units', fontsize=14)
    i+=1
plt.suptitle('Covid-19 Act. Positives - Italy from 25/02/2020', fontsize=16);
plt.show()
