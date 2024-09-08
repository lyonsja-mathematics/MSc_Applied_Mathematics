# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 13:09:08 2016

@author: lyons
"""

from matplotlib.pyplot import *
from math import cos, pi, sqrt
from numpy import *

rc('text',usetex=True)
rc('font',family='serif')

def rho(x,t):
    return 1-cos(x+t)

def rho2(x,t):
    return (1-cos(x*exp(-t)))*exp(-t)

x=linspace(0,2*pi,101)
figure(figsize=(16,6))

subplot(121)
plot(x,rho(x,0),x,rho(x,1),'red',x,rho(x,2),'green',linewidth=3)
legend((r'$t=0$','$t=1$','$t=2$'),loc='best',fontsize=15)
xlabel(r'x',fontsize=20); ylabel(r'\rho',fontsize=20)
title(r'(a)',fontsize=20)
#axis([0,7,0,2.5])

subplot(122)
plot(x,rho2(x,0),x,rho2(x,1),'red',x,rho2(x,2),'green',linewidth=3)
legend(('$t=0$','$t=1$','$t=2$'),loc='best',fontsize=15)
title(r'(b)',fontsize=20)
xlabel(r'x',fontsize=20); ylabel(r'\rho',fontsize=20)

show()