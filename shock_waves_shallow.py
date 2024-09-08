# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 10:23:24 2016

@author: lyons
"""


get_ipython().magic('matplotlib inline')
from numpy import linspace,zeros,exp,arange,sqrt, cos
from warnings import filterwarnings
from math import pi
from matplotlib.pyplot import *

filterwarnings("ignore")

rc('text', usetex=True)
rc('font', family='serif')


def h_0(x):
    value=zeros(npt)
    for i in range(npt):
        if x[i]<0:
            value[i]=1
        if x[i]>=0 and x[i] <= 5:
            value[i]=0.4*x[i]+1
        if x[i] > 5:
            value[i]=1
    return value
'''

def h_0(x):
    return 1-cos(0.5*x)
 
def h_0_s(x):
    if x < 0:
        value=1
    if x>=0 and x<3:
        value=2
    if x >= 3:
        value=-x
    return value
'''
    
def u_01_s(h_star,u_star,T,x):
    return u_star - sqrt(9.81/2)*sqrt( ((h_0_s(x)/h_star) - (h_star/h_0_s(x)))*(h_0_s(x)-h_star))
    
def u_01(h_star,u_star,T,x):
    return u_star + sqrt(9.81/2)*sqrt( ((h_0(x)/h_star) - (h_star/h_0(x)))*(h_0(x)-h_star))
    
def t(h_star,u_star,T,x_star,x):
    return (1/u_01(h_star,u_star,T,x))*(x-x_star)
    
xlim=2*pi
npt=100
dx=xlim/npt

x=linspace(-xlim,xlim,npt)

fig1=figure(figsize=(16,12))



subplot(221)
plot(x,h_0(x),linewidth=3)
axis([-1,6,0.5,3])
xlabel('x',fontsize=20)
ylabel(r'\phi',fontsize=20)
title('(a)',fontsize=20)

subplot(222)
plot(x,u_01(1,0,293,x),linewidth=3)
axis([-1,6,0,6])
xlabel(r'x',fontsize=20)
ylabel(r'u',fontsize=20)
title(r'(b)',fontsize=20)



subplot(223)

for x_star in arange(-3*xlim,3*xlim,dx*15):
    plot(x,t(1,1,293,x_star,x),linewidth=3)
    axis([-1,6,-7,15])
    xlabel(r'x',fontsize=20)
    ylabel(r't',fontsize=20)
    title(r'(c)',fontsize=20)
    
subplot(224)
for x_star in arange(-3*xlim,3*xlim,dx*15):
    plot(x,t(1,1,293,x_star,x),linewidth=3)
    axis([4.85,5.12,-7,15])
    xlabel(r'x',fontsize=20)
    ylabel(r't',fontsize=20)
    title(r'(d)',fontsize=20)
''' 
subplot(132)
for x_star in arange(-xlim,xlim,15*dx):
    plot(x,t(1,1,293,x_star,x),linewidth=3)
    axis([-0.1,0.1,-7,8])
    xlabel('x',fontsize=20)
    ylabel('t',fontsize=20)
    title('(b)',fontsize=20)
    
subplot(133)
for x_star in arange(-xlim,xlim,15*dx):
    plot(x,t(1,1,293,x_star,x),linewidth=3)
    axis([2.9,3.2,-5,10])
    xlabel('x',fontsize=20)
    ylabel('t',fontsize=20)
    title('(c)',fontsize=20)
'''