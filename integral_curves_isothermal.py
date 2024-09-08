# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 22:16:59 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
from numpy import linspace,zeros,exp,arange,sqrt, cos,log
from warnings import filterwarnings
from math import pi
from scipy.optimize import fsolve
from matplotlib.pyplot import *
filterwarnings("ignore")
rc('text',usetex=True)
rc('font', family='serif')

a=sqrt(8.314*293)

def shk2(rho):
    return a*sqrt(((1/0.2)-(1/rho))*(rho-0.2))
    
def rar1(rho):
    return -a*log(rho)
    
rho=linspace(0,2,1000)

plot(rho,rar1(rho),rho,shk2(rho),lw=2)