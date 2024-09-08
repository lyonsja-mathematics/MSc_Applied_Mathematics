# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 16:55:12 2016

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
N=100*(end-start)+1

s=linspace(start,end,N)

def x(s,t):
    x_value=zeros(N)
    for i in range(N):
        x_value[i]=s[i]*t
    return x_value

def rho_r(rho_l,u_l,u_r,a,s):
    return rho_l*exp((u_l-u_r)/a)


def rho(rho_l,u_l,u_r,a,s):
    value=zeros(N)
    for i in range(N):
        if s[i]<=u_l-a:
            value[i]=rho_l
        if s[i] >u_l-a and s[i]<u_r-a:
            value[i]=rho_r(rho_l,u_l,u_r,a,s)*exp((s[i]+a-u_r)/a)+rho_l*exp((u_l-a-s[i])/a)
        if s[i] >=u_r-a:
            value[i]=2*rho_r(rho_l,u_l,u_r,a,s)
    return value


def u(u_l,u_r,a,s):
    value=zeros(N)
    for i in range(N):
        if s[i]<=u_l-a:
            value[i]=2*u_l
        if s[i] > u_l-a and s[i] < u_r-a:
            value[i]=2*s[i]+2*a
        if s[i] >=u_r-a:
            value[i]=2*u_r
    return value

    
figure(figsize=(18,7))

subplot(121)
#a_list=linspace(0,1,11)
#for i in range(len(a_list)):
plot(s,rho(2,-1,1,0.25,s),s,rho(2,-1,1,0.5,s),s,rho(2,-1,1,0.75,s),s,rho(2,-1,1,1,s),linewidth=3)
legend((r'$a=0.25$','$a=0.5$','$a=0.75$','$a=1$'),loc='best',fontsize=20)
title(r'(a)',fontsize=20)
ylabel(r'\rho',fontsize=20)
xlabel(r'x/t',fontsize=20)
axis([-4,4,0,2.1])

subplot(122)
a_list=linspace(0,1,5)
#for i in range(len(a_list)):
plot(s,u(-1,1,0.25,s),s,u(-1,1,0.50,s),s,u(-1,1,0.75,s),s,u(-1,1,1,s),linewidth=3)
legend((r'$a=0.25$','$a=0.5$','$a=0.75$','$a=1$'),loc='best',fontsize=20)
axis([-4,4,-2.1,2.1])
ylabel(r'u',fontsize=20)
xlabel(r'x/t',fontsize=20)
title(r'(b)',fontsize=20)

show()