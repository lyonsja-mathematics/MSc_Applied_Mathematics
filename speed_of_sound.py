# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 19:06:44 2016

@author: lyons
"""

from matplotlib.pyplot import *
from numpy import linspace, sqrt,cos,sin,log
import warnings
warnings.filterwarnings("ignore")

rc('text', usetex=True)
rc('font', family='serif')

def c(r,kappa):
    gamma=1.31
    return sqrt(gamma*kappa)*r**((gamma-1)/2)/340

rlim=4
r=linspace(0,rlim,101)
figure(figsize=(7,5))
kappa=10000
while kappa <= 100000:
    plot(r,c(r,kappa),linewidth=3)
    kappa = kappa+10000
xlabel(r'\rho',fontsize=20); ylabel(r'c',fontsize=20)
show()