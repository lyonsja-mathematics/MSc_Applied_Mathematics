# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 18:38:10 2016

@author: lyons
"""

from matplotlib.pyplot import *
from numpy import linspace, sqrt
#from math import sqrt
import warnings
warnings.filterwarnings("ignore")
#from math import log

rc('text', usetex=True)
rc('font', family='serif')

#u0=2
phi0=2

def v1(phi,u0):
    return u0 +2*(sqrt(phi0)-sqrt(phi))

def v2(phi,u0):
    return u0 +2*(sqrt(phi)-sqrt(phi0))


philim=4
phi=linspace(0,philim, 101)

figure(figsize=(16,6))

subplot(121)
for u0 in range(10):
    plot(phi,v1(phi,u0),linewidth=3)
xlabel(r'\phi',fontsize=20)
ylabel(r'u',fontsize=20)
axis([0,4,0,10])
title(r'(a)',fontsize=20)

subplot(122)
for u0 in range(10):
    plot(phi,v2(phi,u0),linewidth=3)
xlabel(r'\phi',fontsize=20)
ylabel(r'u',fontsize=20)
title(r'(b)',fontsize=20)
axis([0,4,0,8])
show()