# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 18:51:25 2016

@author: lyons
"""

from matplotlib.pyplot import *
from numpy import linspace, sqrt,cos,sin
#from math import sqrt
import warnings
warnings.filterwarnings("ignore")
#from math import log

rc('text', usetex=True)
rc('font', family='serif')

c=0.8

def rho(x,t):
    return 1 - cos(x-t)*cos(c*t)+(1/c)*sin(x-t)*sin(c*t)

def u(x,t):
    return 1-c*sin(x-t)*sin(c*t)+cos(x-t)*cos(c*t)

x=linspace(0,7)

figure(figsize=(16,6))

subplot(121)
#for t in range(3):
plot(x,rho(x,0),x,rho(x,1),x,rho(x,2),linewidth=3)
xlabel(r'x',fontsize=20); ylabel(r'\rho',fontsize=20)
legend(('$t=0$','$t=1$','$t=2$'),loc='best',fontsize=15)

subplot(122)
#for t in range(3):
plot(x,u(x,0),x,u(x,1),x,u(x,2),linewidth=3)
xlabel(r'x',fontsize=20); ylabel(r'u',fontsize=20)
legend(('$t=0$','$t=1$','$t=2$'),loc='best',fontsize=15)

show()