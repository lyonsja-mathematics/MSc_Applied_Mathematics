# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 18:05:18 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import warnings
warnings.filterwarnings("ignore")

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

npt=1000
a=50
xmin=-1
xmax=1

x=np.linspace(xmin,xmax,npt)

def H(x):
    v=np.zeros(npt)
    for i in range(npt):
        if x[i]<0:
            v[i]=0
        if x[i]>0:
            v[i]=1
    return v
    
def rho_0(rhol,rhor):
    return (rhol+rhor)/2

def u_0(ul,ur):
    return (ul+ur)/2

def w_1_l(rhol,ul,rhor,ur):
    rho0=rho_0(rhol,rhor)
    return (ul/2) - (a*rhol)/(2*rho0)
    
def w_2_l(rhol,ul,rhor,ur):
    rho0=rho_0(rhol,rhor)
    return (ul/2) + (a*rhol)/(2*rho0)
    
def w_1_r(rhol,ul,rhor,ur):
    rho0=rho_0(rhol,rhor)
    return (ur/2) - (a*rhor)/(2*rho0)
    
def w_2_r(rhol,ul,rhor,ur):
    rho0=rho_0(rhol,rhor)
    return (ur/2) + (a*rhor)/(2*rho0)
    
def u_lin(rhol,ul,rhor,ur,t):
    
    u0=u_0(ul,ur)
    rho0=rho_0(rhol,rhor)
    
    lambda1=u0-a
    lambda2=u0+a
    
    w1r=w_1_r(rhol,ul,rhor,ur)
    w2r=w_2_r(rhol,ul,rhor,ur)
    
    w1l=w_1_l(rhol,ul,rhor,ur)
    w2l=w_2_l(rhol,ul,rhor,ur)
    
    ta = t*np.ones(npt)
    
    return ul + H(x-lambda1*ta)*(w1r-w1l)+H(x-lambda2*ta)*(w2r-w2l)
'''    
def plot_lin_sol_u(rhol,ul,rhor,ur,t):
    
    
    plt.plot(x,u(rhol,ul,rhor,ur,t),lw=3)

    umax=max(u(rhol,ul,rhor,ur,t))
    umin=min(u(rhol,ul,rhor,ur,t))
    udiff=umax-umin
    plt.axis([xmin,xmax,umin-0.1*udiff,umax+0.1*udiff])
    plt.xlabel(r'$x$',fontsize=20)
    plt.ylabel(r'$u$',fontsize=20)
    plt.title(r'$t=%.3f$' % t,fontsize=20)
'''  
def rho_lin(rhol,ul,rhor,ur,t):
    
    u0=u_0(ul,ur)
    rho0=rho_0(rhol,rhor)
    
    lambda1=u0-a
    lambda2=u0+a
    
    w1r=w_1_r(rhol,ul,rhor,ur)
    w2r=w_2_r(rhol,ul,rhor,ur)
    
    w1l=w_1_l(rhol,ul,rhor,ur)
    w2l=w_2_l(rhol,ul,rhor,ur)
    
    ta = t*np.ones(npt)
    
    return rhol - H(x-lambda1*ta)*(w1r-w1l)*(rho0/a) + H(x-lambda2*ta)*(w2r-w2l)*(rho0/a)
    
def plot_lin_sol_rho(rhol,ul,rhor,ur,t):
    
    
    plt.plot(x,rho(rhol,ul,rhor,ur,t),lw=3)

    rhomax=max(rho(rhol,ul,rhor,ur,t))
    rhomin=min(rho(rhol,ul,rhor,ur,t))
    rhodiff=rhomax-rhomin
    plt.axis([xmin,xmax,rhomin-0.1*rhodiff,rhomax+0.1*rhodiff])
    plt.xlabel(r'$x$',fontsize=20)
    plt.ylabel(r'$\rho$',fontsize=20)
    plt.title(r'$t=%.3f$' % t,fontsize=20)


def plot_rho_u(rhol,ul,rhor,ur,t):
    
    plt.figure(figsize=(15,6))

    plt.subplot(121)
    plot_lin_sol_rho(rhol,ul,rhor,ur,t)

    plt.subplot(122)
    plot_lin_sol_u(rhol,ul,rhor,ur,t)

#fot t in [0,0.01,]
plot_rho_u(1,0,2,0,0.01)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    