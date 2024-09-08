# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 23:15:32 2016

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

a=np.sqrt(8.314*293)
npt=1000

def rho_0(rhol,ul,rhor,ur):
    return (ul-ur)/((a/rhor)+(a/rhol))
    
def u_0(rhol,ul,rhor,ur):
    return (rhol*ul +rhor*ur)/(rhor+rhol)
    
def rho_init(rhol,rhor,x):
    value=np.zeros(npt)
    for i in range(npt):
        if x[i]<=0:
            value[i]=rhol
        if x[i]>0:
            value[i]=rhor
    return value
    
def u_init(ul,ur,x):
    value=np.zeros(npt)
    for i in range(npt):
        if x[i]<=0:
            value[i]=ul
        if x[i]>0:
            value[i]=ur
    return value
        
def w1(rhol,ul,rhor,ur,x,t):
    
    rho0=rho_0(rhol,ul,rhor,ur)
    lambda1=u_0(rhol,ul,rhor,ur)-a

    return (1/(2*rho0))*(rho0*u_init(ul,ur,x-lambda1*t) - a*rho_init(rhol,rhor,x-lambda1*t))
    
    
def w2(rhol,ul,rhor,ur,x,t):
    
    rho0=rho_0(rhol,ul,rhor,ur)
    lambda2=u_0(rhol,ul,rhor,ur)+a

    return (1/(2*rho0))*(rho0*u_init(ul,ur,x-lambda2*t) + a*rho_init(rhol,rhor,x-lambda2*t))

    
def rho_lin(rhol,ul,rhor,ur,x,t):
    
    rho0=rho_0(rhol,ul,rhor,ur)
    return (rho0/a)*(w2(rhol,ul,rhor,ur,x,t)-w1(rhol,ul,rhor,ur,x,t))
    
def u_lin(rhol,ul,rhor,ur,x,t):
    return w2(rhol,ul,rhor,ur,x,t) + w1(rhol,ul,rhor,ur,x,t)
    
def plot_lin_sol_u(rhol,ul,rhor,ur,t):
    
    xmin=-5
    xmax=5
    npt=1000
    x=np.linspace(xmin,xmax,npt)

    #

    #plt.subplot(121)
    plt.plot(x,u(rhol,ul,rhor,ur,x,t),lw=3)

    umax=max(u(rhol,ul,rhor,ur,x,t))
    umin=min(u(rhol,ul,rhor,ur,x,t))
    udiff=umax-umin
    plt.axis([xmin,xmax,umin-0.1*udiff,umax+0.1*udiff])
    plt.xlabel(r'$x$',fontsize=20)
    plt.ylabel(r'$u$',fontsize=20)
    plt.title(r'$t=%.3f$' % t,fontsize=20)

#plot_lin_sol_u(1,-1,0.1,1,x,0.05)


def plot_lin_sol_rho(rhol,ul,rhor,ur,t):
    
    xmin=-5
    xmax=5
    npt=1000
    x=np.linspace(xmin,xmax,npt)

    #plt.figure(figsize=(15,6))

    #plt.subplot(121)
    plt.plot(x,rho(rhol,ul,rhor,ur,x,t),lw=3)

    rhomax=max(rho(rhol,ul,rhor,ur,x,t))
    rhomin=min(rho(rhol,ul,rhor,ur,x,t))
    rhodiff=rhomax-rhomin
    plt.axis([xmin,xmax,rhomin-0.1*rhodiff,rhomax+0.1*rhodiff])
    plt.xlabel(r'$x$',fontsize=20)
    plt.ylabel(r'$\rho$',fontsize=20)
    plt.title(r'$t=%.3f$' % t,fontsize=20)

#plot_lin_sol_rho(1,-1,2,1,x,0.02)
def plot_rho_u(rhol,ul,rhor,ur,t):
    
    plt.figure(figsize=(15,6))

    plt.subplot(121)
    plot_lin_sol_rho(rhol,ul,rhor,ur,t)

    plt.subplot(122)
    plot_lin_sol_u(rhol,ul,rhor,ur,t)
    

    
    
plot_rho_u(1,-0.1,0.1,0.1,0.01)
    
