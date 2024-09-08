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
rc('text', usetex=True)
rc('font', family='serif')


def rho_0(x):
    value=zeros(npt)
    for i in range(npt):
        if x[i]<0:
            value[i]=1
        #if x[i]>=0 and x[i] < 3:
         #   value[i]=2
        if x[i] >= 0:
            value[i]=2
    return value
'''    
def rho_1(x):
    value=zeros(npt)
    for i in range(npt):
        if x[i]<0:
            value[i]=1.273
        #if x[i]>=0 and x[i] < 3:
         #   value[i]=2
        if x[i] >= 0:
            value[i]=1
    return value    

   
def rho_0_s(x):
    if x < 0:
        value=1
    #if x>=0 and x<3:
     #   value=2
    if x >= 0:
        value=2
    return value
    
def u_01_s(rho_star,u_star,T,x):
    return u_star - sqrt(8.314*T)*sqrt( ((1/rho_star) - (1/rho_0_s(x)))*(rho_0_s(x)-rho_star))
'''   
def u_01(rho_star,u_star,T,x):
    return u_star - sqrt(8.314*T)*sqrt( ((1/rho_star) - (1/rho_0(x)))*(rho_0(x)-rho_star))

def t(rho_star,u_star,T,x_star,x):
    return (1/u_01(rho_star,u_star,T,x))*(x-x_star)  
    
'''   
def u_02(rho_star,u_star,T,x):
    return u_star - sqrt(8.314*T)*sqrt( ((1/rho_star) - (1/rho_1(x)))*(rho_1(x)-rho_star))
    

    
def t_1(rho_star,u_star,T,x_star,x):
    return (1/u_02(rho_star,u_star,T,x))*(x-x_star)
'''    
xlim=2*pi
npt=100
dx=xlim/npt

x=linspace(-xlim,xlim,npt)

#fig1=figure(figsize=(6,6))


'''
subplot(121)
plot(x,rho_0(x),linewidth=3)
axis([-6,6,0.8,2.1])
xlabel('x',fontsize=20)
ylabel('Density',fontsize=20)
title('(a)',fontsize=20)

subplot(122)
plot(x,u_01(1,1,293,x),linewidth=3)
axis([-6,6,-35,5])
xlabel('x',fontsize=20)
ylabel('u',fontsize=20)
title('(b)',fontsize=20)
'''

fig1=figure(figsize=(8,6))
#subplot(131)

#subplot(121)
for x_star in arange(-4*xlim,4*xlim,dx*20):
    plot(x,t(1.1,0,293,x_star,x),linewidth=3)
    axis([-6,6,-0.5,0.5])
    xlabel('$x$',fontsize=30)
    ylabel(r'$t$',fontsize=30)
    #title('(a)',fontsize=20)
'''     
subplot(122)
for x_star in arange(-4*xlim,4*xlim,dx*20):
    plot(x,t(1.75,0,293,x_star,x),linewidth=3)
    axis([-0.1,0.1,-0.5,0.5])
    xlabel('x',fontsize=20)
    ylabel('t',fontsize=20)
    title('(b)',fontsize=20)
   
subplot(122)
for x_star in arange(-4*xlim,4*xlim,dx*20):
    plot(x,t_1(1,0,293,x_star,x),linewidth=3)
    #axis([-6,6,-0.4,0.4])
    xlabel('x',fontsize=20)
    ylabel('t',fontsize=20)
'''
'''   
subplot(132)
for x_star in arange(-3*xlim,3*xlim,15*dx):
    plot(x,t(1,1,293,x_star,x),linewidth=3)
    axis([-0.1,0.1,-7,8])
    xlabel('x',fontsize=20)
    ylabel('t',fontsize=20)
    title('(b)',fontsize=20)
    
subplot(133)
for x_star in arange(-3*xlim,3*xlim,10*dx):
    plot(x,t(1,1,293,x_star,x),linewidth=3)
    axis([0.059,0.064,-0.2,0.2])
    xlabel('x',fontsize=20)
    ylabel('t',fontsize=20)
    title('(c)',fontsize=20)
'''
show()