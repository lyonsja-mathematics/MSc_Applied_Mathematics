# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 12:28:09 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
from numpy import linspace,zeros,exp
from warnings import filterwarnings
from matplotlib.pyplot import *

filterwarnings("ignore")

rc('text', usetex=True)
rc('font', family='serif')


start=-10
end=10
N=10*(end-start)+1

s=linspace(start,end,N)

def x(s,t):
    x_value=zeros(N)
    for i in range(N):
        x_value[i]=s[i]*t
    return x_value


def rho_1(s):
    value=zeros(N)
    for i in range(N):
        if s[i] <= -1.7:
            value[i]=2
        if s[i] > -1.7 and s[i] < 0.3:
            value[i]=2*exp((-17-10*s[i])/7)
        if s[i] >= 0.3:
            value[i]=0.115
    return value

def u(s):
    u_value=zeros(N)
    for i in range(N):
        if s[i] <= -1.7:
            u_value[i]=-1
        if s[i] > -1.7 and s[i] < 0.3:
            u_value[i]=s[i] +0.7
        if s[i] >= 0.3:
            u_value[i]=1
    return u_value
    
    
def rho_2(s):
    value=zeros(N)
    for i in range(N):
        if s[i] <= -1.7:
            value[i]=2
        if s[i] > -1.7 and s[i] < 0.3:
            value[i]=34.82*exp((-3+10*s[i])/7)
        if s[i] >= 0.3:
            value[i]=34.82
    return value
    



figure(figsize=(20,16))

subplot(221)
#for t in range(1,15):
plot(x(s,1),rho_1(s)+rho_2(s),x(s,1),rho_1(s),x(s,1),rho_2(s),linewidth=3)
xlabel(r'x/t',fontsize=20)
ylabel(r'\rho',fontsize=20)
#axis([-10,10,0,2.1])
legend((r'Superposition','Field-1','Field-2',),loc='best',fontsize=15)
title(r'(a)',fontsize=20)

subplot(222)
#for t in range(1,15):
plot(x(s,1),rho_1(s)+rho_2(s),x(s,1),rho_1(s),x(s,1),rho_2(s),linewidth=3)
xlabel(r'x/t',fontsize=20)
ylabel(r'\rho',fontsize=20)
axis([-10,10,34.6,35])
legend((r'Superposition','Field-1','Field-2',),loc='best',fontsize=15)
title(r'(b)',fontsize=20)

subplot(223)
#for t in range(1,15):
plot(x(s,1),u(s)+u(s),x(s,1),u(s),x(s,1),u(s),linewidth=3)
xlabel(r'x/t',fontsize=20)
ylabel(r'u',fontsize=20)
axis([-10,10,-2.1,2.1])
legend((r'Superposition','Field-1','Field-2',),loc='best',fontsize=15)
title(r'(c)',fontsize=20)