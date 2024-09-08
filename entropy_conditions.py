# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 09:24:59 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
from numpy import linspace,zeros,sqrt
from warnings import filterwarnings
from matplotlib.pyplot import *

filterwarnings("ignore")

rc('text', usetex=True)
rc('font', family='serif')

rho=linspace(0,5,51)


def rho_0(x):
    value=zeros(51)
    for i in range(51):
        if x[i]<0:
            value[i]=1
        #if x[i]>=0 and x[i] < 3:
         #   value[i]=2
        if x[i] >= 0:
            value[i]=2
    return value

def u_1(rho_star,u_star,kappa,gamma,rho):
    return u_star -sqrt(kappa)*sqrt(((1/rho_star)-(1/rho))*(rho**gamma-rho_star**gamma))
    
def u_2(rho_star,u_star,kappa,gamma,rho):
    return u_star +sqrt(kappa)*sqrt(((1/rho_star)-(1/rho))*(rho**gamma-rho_star**gamma))


gamma_list=[1.67,1.4,1.33,1]
T_list=[253,273,293,373]
kappa_list = [1,2,3]

figure(figsize=(16,16))

subplot(221)
for i in range(len(gamma_list)):
    plot(rho,u_1(1,0,1,gamma_list[i],rho),linewidth=3)
xlabel(r'\rho',fontsize=20)
ylabel(r'u',fontsize=20)
axis([1,4,-6,0])
legend((r'$\gamma = 1.67$','$\gamma = 1.4$','$\gamma = 1.33$','$\gamma = 1$'), fontsize=20, loc='best')
title(r'(a)',fontsize=20)

subplot(222)
for i in range(len(gamma_list)):
    plot(rho,u_2(1,0,1,gamma_list[i],rho),linewidth=3)
xlabel(r'\rho',fontsize=20)
ylabel(r'u',fontsize=20)
axis([1,4,0,6])
legend((r'$\gamma = 1.67$','$\gamma = 1.4$','$\gamma = 1.33$','$\gamma = 1$'), fontsize=20, loc='best')
title(r'(b)',fontsize=20)

subplot(223)
for i in range(len(kappa_list)):
    plot(rho,u_1(1,0,kappa_list[i],1.4,rho),linewidth=3)
xlabel(r'\rho',fontsize=20)
ylabel(r'u',fontsize=20)
axis([1,4,-6,0])
legend((r'$\kappa = 1$','$\kappa = 2$','$\kappa = 3$'), fontsize=20, loc='best')
title(r'(c)',fontsize=20)

subplot(224)
for i in range(len(kappa_list)):
    plot(rho,u_2(1,0,kappa_list[i],1.4,rho),linewidth=3)
xlabel(r'\rho',fontsize=20)
ylabel(r'u',fontsize=20)
axis([1,4,0,6])
legend((r'$\kappa = 1$','$\kappa = 2$','$\kappa = 3$'), fontsize=20, loc='best')
title(r'(d)',fontsize=20)





