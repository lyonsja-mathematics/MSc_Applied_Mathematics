# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 20:44:20 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
from numpy import linspace,zeros,sqrt,arange
from warnings import filterwarnings
from matplotlib.pyplot import *

filterwarnings("ignore")

rc('text', usetex=True)
rc('font', family='serif')

def rho(x):
    value=zeros(npt)
    for i in range(npt):
        if x[i]<0:
            value[i]=2
        if x[i] >= 0:
            value[i]=1
    return value

def u_1(rho_star,u_star,kappa,gamma,x):
    return u_star - (2/gamma-1)*sqrt(kappa*gamma)*(rho(x)**((gamma-1)/2)-rho_star**((gamma-1)/2))
    
def u_2(rho_star,u_star,kappa,gamma,x):
    return u_star + (2/gamma-1)*sqrt(kappa*gamma)*(rho(x)**((gamma-1)/2)-rho_star**((gamma-1)/2))
    
def t_1(rho_star,u_star,kappa,gamma,x_star,x):
    return (1/u_1(rho_star,u_star,kappa,gamma,x))*(x-x_star)
    
def t_2(rho_star,u_star,kappa,gamma,x_star,x):
    return (1/u_2(rho_star,u_star,kappa,gamma,x))*(x-x_star)
    
xlim=6
npt=1000
dx=xlim/npt

x=linspace(-xlim,xlim,npt)

fig1=figure(figsize=(9,9))

#subplot(121)
for x_star in arange(-4*xlim,4*xlim,dx*200):
    plot(x,t_1(0.8,0,1,1.4,x_star,x),linewidth=3)
    axis([-6,6,-100,100])
    xlabel('x',fontsize=20)
    ylabel('t',fontsize=20)
    #title('(a)',fontsize=20)
