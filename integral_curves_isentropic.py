# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 18:57:43 2016

@author: lyons
"""

from matplotlib.pyplot import *
from numpy import linspace, sqrt,cos,sin,log
#from math import sqrt
import warnings
warnings.filterwarnings("ignore")
#from math import log

rc('text', usetex=True)
rc('font', family='serif')


c=0.8; u0=1; r0=3

def lv1(r):
    return u0 +c - (c/r0)*r

def lv2(r):
    return u0 -c + (c/r0)*r

def v1(r):
    return u0 - c*log(r/r0)

def v2(r):
    return u0 + c*log(r/r0)

rholim=6
rho=linspace(0,rholim, 101)

figure(figsize=(16,6))

subplot(121)
plot(rho,lv1(rho), rho, v1(rho),linewidth=3)
xlabel(r'\rho',fontsize=20); ylabel(r'u',fontsize=20)
legend(('Linearised model', 'Non linearised'), loc='best',fontsize=15)
title('(a)',fontsize=20)

subplot(122)
plot(rho,lv2(rho), rho, v2(rho),linewidth=3)
xlabel(r'\rho',fontsize=20); ylabel(r'u',fontsize=20)
legend(('Linearised system', 'Non linearised'), loc='best',fontsize=15)
title('(b)',fontsize=20)

