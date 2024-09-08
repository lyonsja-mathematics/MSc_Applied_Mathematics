# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 16:31:10 2016

@author: lyons
"""


get_ipython().magic('matplotlib inline')
from numpy import linspace,zeros,exp,arange,sqrt
from warnings import filterwarnings
from math import pi, cos
from matplotlib.pyplot import *

filterwarnings("ignore")

def rho_0(x):
    value=zeros(npt)
    for i in range(npt):
        if x[i]<0:
            value[i]=1
        if x[i] >= 0:
            value[i]=2
    return value
    

def u_1(rho_star,u_star,T,x):
    return u_star - sqrt(8.314*T)*sqrt( ((1/rho_star) - (1/rho_0(x)))*(rho_0(x)-rho_star))
    
def u_2(rho_star,u_star,T,x):
    return u_star + sqrt(8.314*T)*sqrt( ((1/rho_star) - (1/rho_0(x)))*(rho_0(x)-rho_star))
    
def t_1(rho_star,u_star,T,x_star,x):
    return (1/u_1(rho_star,u_star,T,x))*(x-x_star)
    
def t_2(rho_star,u_star,T,x_star,x):
    return (1/u_2(rho_star,u_star,T,x))*(x-x_star)
    
xlim=2*pi
npt=1000
dx=xlim/npt

x=linspace(-xlim,xlim,npt)

fig1=figure(figsize=(18,9))


subplot(121)
for x_star in arange(-5*xlim,5*xlim,dx*250):
    plot(x,t_1(0.5,0,293,x_star,x),linewidth=3)
    axis([-6,6,-0.3,0.3])
    xlabel('x',fontsize=20)
    ylabel('t',fontsize=20)
    title('(a)',fontsize=20)
    
subplot(122)
for x_star in arange(-5*xlim,5*xlim,dx*250):
    plot(x,t_2(0.5,0,293,x_star,x),linewidth=3)
    axis([-6,6,-0.3,0.3])
    xlabel('x',fontsize=20)
    ylabel('t',fontsize=20)
    title('(b)',fontsize=20)

show()