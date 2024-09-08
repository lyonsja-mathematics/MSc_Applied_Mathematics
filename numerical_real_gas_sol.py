# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 22:57:58 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math as m
import warnings
warnings.filterwarnings("ignore")

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


h=0.05
r=0.001
k=r*h


gamma=1.4
a=0.1382

n=int(2/h)
m=int(1/k)

x=np.linspace(-1,1,n)

rho=np.ones((n,m))
for i in range(n):
    if i <= n/2:
        rho[i][0]=3
    if i > n/2:
        rho[i][0]=1


p=np.ones((n,m))
for i in range(n):
    if i <= n/2:
        p[i][0]=3
    if i > n/2:
        p[i][0]=1

u=np.ones((n,m))
for i in range(n):
    u[i][0]=0

c2=np.zeros((n,m))
for i in range(n):
    for j in range(m):
        c2[i][j]=gamma*p[i][j]/rho[i][j]

for j in range(m-1):
    for i in range(n-1):
            
        p[i][j+1] = p[i][j] + r*(gamma-1)*c2[i][j]*u[i][j]*(rho[i+1][j]-rho[i][j])-r*gamma*p[i][j]*(u[i+1][j]-u[i][j]) + r*a*(gamma-1)*rho[i][j]*rho[i][j]*(u[i+1][j]-u[i][j]) - r*gamma*u[i][j]*(p[i+1][j]-p[i][j])
        
        rho[i][j+1]=rho[i][j]-r*rho[i][j]*(u[i+1][j]-u[i][j])-r*u[i][j]*(rho[i+1][j]-rho[i][j])
        
        u[i][j+1]=u[i][j]-r*u[i][j]*(u[i+1][j]-u[i][j])-(r/rho[i][j])*(p[i+1][j]-p[i][j])


        
        
def pi(t):
    j=int(t/k)
    val=np.zeros(n)
    for i in range(n):
        val[i]=p[i][j]
    return val
 
def plot_p(t):
    
    plt.plot(x,pi(t),lw=2)
    
    fts=30
    pmax=max(pi(t))
    pmin=min(pi(t))
    pdiff=pmax-pmin
    plt.axis([-1,1,pmin-0.1*pdiff,pmax+0.1*pdiff])
    plt.xlabel(r'$x$',fontsize=fts)
    plt.ylabel(r'$p$',fontsize=fts)
    
def rhoi(t):
    j=int(t/k)
    val=np.zeros(n)
    for i in range(n):
        val[i]=rho[i][j]
    return val
 
def plot_rho(t):
    
    plt.plot(x,rhoi(t),lw=2)
    
    fts=30
    pmax=max(rhoi(t))
    pmin=min(rhoi(t))
    pdiff=pmax-pmin
    plt.axis([-1,1,0,3.2])
    plt.xlabel(r'$x$',fontsize=fts)
    plt.ylabel(r'$\rho$',fontsize=fts)
    
def ui(t):
    j=int(t/k)
    val=np.zeros(n)
    for i in range(n):
        val[i]=u[i][j]
    return val
 
def plot_u(t):
    
    plt.plot(x,ui(t),lw=2)
    
    fts=30
    pmax=max(ui(t))
    pmin=min(ui(t))
    pdiff=pmax-pmin
    plt.axis([-1,1,pmin-0.1*pdiff,pmax+0.1*pdiff])
    plt.xlabel(r'$x$',fontsize=fts)
    plt.ylabel(r'$u$',fontsize=fts)
    
plot_rho(0.5)
#print(rho)