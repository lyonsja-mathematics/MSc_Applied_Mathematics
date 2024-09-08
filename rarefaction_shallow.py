# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 08:20:08 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
from numpy import linspace,zeros,sqrt
from warnings import filterwarnings
from matplotlib.pyplot import *

filterwarnings("ignore")

rc('text', usetex=True)
rc('font', family='serif')

start=-10
end=10
N=100*(end-start)+1

s=linspace(start,end,N)

def x(s,t):
    x_value=zeros(N)
    for i in range(N):
        x_value[i]=s[i]*t
    return x_value

def phi_2(s):
    value=zeros(N)
    for i in range(N):
        if s[i] <= -0.45:
            value[i]=0.2
        if s[i] > -0.45 and s[i] < 2.55:
            value[i]=(sqrt(2.1)+(1/3)*(s[i]-2.55))**2
        if s[i] >= 2.55:
            value[i]=2.1
    return value

def phi_1(s):
    value=zeros(N)
    for i in range(N):
        if s[i] <= -0.45:
            value[i]=6
        if s[i] > -0.45 and s[i] < 2.55:
            value[i]=(sqrt(6)-(1/3)*(s[i]+0.45))**2
        if s[i] >= 2.55:
            value[i]=2.1
    return value


def u(s):
    u_value=zeros(N)
    for i in range(N):
        if s[i] <= -0.45:
            u_value[i]=-1
        if s[i] > -0.45 and s[i] < 2.55:
            u_value[i]=(2/3)*s[i] -0.7
        if s[i] >=2.55:
            u_value[i]=1
    return u_value


figure(figsize=(16,6))

subplot(121)
#for t in range(1,15):
plot(x(s,1),phi_1(s),x(s,1),phi_2(s),linewidth=3)
legend((r'Field-1','Field-2'),loc='best')
xlabel(r'x',fontsize=20)
ylabel(r'\phi',fontsize=20)
axis([start,end,0,6.1])
title(r'(a)',fontsize=20)

subplot(122)
#for t in range(1,15):
plot(x(s,1),u(s),linewidth=3)
#legend(('2','1','0.5','0','-0.5'),loc='best')
xlabel(r'x',fontsize=20)
ylabel(r'u',fontsize=20)
axis([start,end,-1.1,1.1])
title(r'(b)',fontsize=20)

show()