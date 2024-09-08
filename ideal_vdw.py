# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 14:50:28 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
from numpy import linspace,zeros,exp,arange
from warnings import filterwarnings
from math import pi, cos
from matplotlib.pyplot import *

filterwarnings("ignore")

rc('text',usetex=True)
rc('font',family='serif')

def ideal(rho,T):
    return rho*8.314*T
    
def vdw(rho,T,a,b):
    return (rho*8.314*T)/(1-rho*b) - a*rho**2

figure(figsize=(10,10))

rho=linspace(0,5,501)

T=293

plot(rho,ideal(rho,T),rho,vdw(rho,T,3.66,0.043),lw=3)
xlabel(r'\rho',fontsize=20)
ylabel(r'p',fontsize=20)
legend((r'Ideal','van der Waals'),fontsize=15,loc='best')