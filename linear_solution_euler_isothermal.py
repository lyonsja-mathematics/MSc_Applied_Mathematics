# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 10:39:07 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import warnings
import math as m
warnings.filterwarnings("ignore")

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

npt=1000
a=50
xmin=-1
xmax=1

fts=30

x=np.linspace(xmin,xmax,npt)

a=m.sqrt(8.314*293)

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

def p_0(pl,pr):
    return (pl+pr)/2


def c(rhol,cl,rho):
    gamma=1.4
    return cl*(rho/rhol)**((gamma-1)/2)
    

def w_l(rhol,ul,pl,cl,rhor,ur,pr):
    rho0=rho_0(rhol,rhor)
    c0=a
    
    return rhol-(pl/(c0**2)), 0.5*ul - (pl/(2*rho0*c0)), 0.5*ul + (pl/(2*rho0*c0)), rhol-(pl/(a**2)), 0.5*ul - (pl/(2*rho0*a)), 0.5*ul + (pl/(2*rho0*a))
    
            
            
def w_r(rhol,ul,pl,cl,rhor,ur,pr):
    rho0=rho_0(rhol,rhor)
    c0=a
    
    return rhor-(pr/(c0**2)), 0.5*ur - (pr/(2*rho0*c0)), 0.5*ur + (pr/(2*rho0*c0)), rhor-(pr/(a**2)), 0.5*ur - (pr/(2*rho0*a)), 0.5*ur + (pr/(2*rho0*a))
   
def u_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t):
    
    u0=u_0(ul,ur)
    rho0=rho_0(rhol,rhor)
    c0=a
    
    lambda1=u0
    lambda2=u0-c0
    lambda3=u0+c0
    
    w1r=w_r(rhol,ul,pl,cl,rhor,ur,pr)[0]
    w2r=w_r(rhol,ul,pl,cl,rhor,ur,pr)[1]
    w3r=w_r(rhol,ul,pl,cl,rhor,ur,pr)[2]
    
    w1l=w_l(rhol,ul,pl,cl,rhor,ur,pr)[0]
    w2l=w_l(rhol,ul,pl,cl,rhor,ur,pr)[1]
    w3l=w_l(rhol,ul,pl,cl,rhor,ur,pr)[2]
    
    ta = t*np.ones(npt)
    
    return ul +H(x-lambda2*ta)*(w2r-w2l) + H(x-lambda3*ta)*(w3r-w3l)
    
def plot_u_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t):
    
    
    plt.plot(x,u_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t),lw=3)

    umax=max(u_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t))
    umin=min(u_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t))
    udiff=umax-umin
    
    plt.axis([xmin,xmax,umin-0.1*udiff,umax+0.1*udiff])
    plt.xlabel(r'$x$',fontsize=fts)
    plt.ylabel(r'$u$',fontsize=fts)
    plt.title(r'$t=%.3f$' % t,fontsize=fts)
    
#plot_u_lin(1,0,1,1,0.1,0,2,0.03)

  
def rho_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t):
    
    u0=u_0(ul,ur)
    rho0=rho_0(rhol,rhor)
    c0=a
    
    lambda1=u0
    lambda2=u0-c0
    lambda3=u0+c0
    
    w1r=w_r(rhol,ul,pl,cl,rhor,ur,pr)[0]
    w2r=w_r(rhol,ul,pl,cl,rhor,ur,pr)[1]
    w3r=w_r(rhol,ul,pl,cl,rhor,ur,pr)[2]
    
    w1l=w_l(rhol,ul,pl,cl,rhor,ur,pr)[0]
    w2l=w_l(rhol,ul,pl,cl,rhor,ur,pr)[1]
    w3l=w_l(rhol,ul,pl,cl,rhor,ur,pr)[2]
    
    ta = t*np.ones(npt)
    
    return rhol +H(x-lambda1*ta)*(w1r-w1l) - H(x-lambda2*ta)*(w2r-w2l)*(rho0/c0) + H(x-lambda3*ta)*(w3r-w3l)*(rho0/c0)
    
def plot_rho_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t):
    
    
    plt.plot(x,rho_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t),lw=3)

    rhomax=max(rho_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t))
    rhomin=min(rho_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t))
    rhodiff=rhomax-rhomin
    
    plt.axis([xmin,xmax,rhomin-0.1*rhodiff,rhomax+0.1*rhodiff])
    plt.xlabel(r'$x$',fontsize=fts)
    plt.ylabel(r'$\rho$',fontsize=fts)
    plt.title(r'$t=%.3f$' % t,fontsize=fts)
    
def p_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t):
    
    u0=u_0(ul,ur)
    rho0=rho_0(rhol,rhor)
    c0=a
    
    lambda1=u0
    lambda2=u0-c0
    lambda3=u0+c0
    
    w1r=w_r(rhol,ul,pl,cl,rhor,ur,pr)[0]
    w2r=w_r(rhol,ul,pl,cl,rhor,ur,pr)[1]
    w3r=w_r(rhol,ul,pl,cl,rhor,ur,pr)[2]
    
    w1l=w_l(rhol,ul,pl,cl,rhor,ur,pr)[0]
    w2l=w_l(rhol,ul,pl,cl,rhor,ur,pr)[1]
    w3l=w_l(rhol,ul,pl,cl,rhor,ur,pr)[2]
    
    ta = t*np.ones(npt)
    
    return pl - H(x-lambda2*ta)*(w2r-w2l)*rho0*c0 + H(x-lambda3*ta)*(w3r-w3l)*rho0*c0

def plot_p_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t):
    
    
    plt.plot(x,p_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t),lw=3)

    pmax=max(p_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t))
    pmin=min(p_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t))
    pdiff=pmax-pmin
    fts=30
    plt.axis([xmin,xmax,pmin-0.1*pdiff,pmax+0.1*pdiff])
    plt.xlabel(r'$x$',fontsize=fts)
    plt.ylabel(r'$p$',fontsize=fts)
    plt.title(r'$t=%.3f$' % t,fontsize=fts)
    
def E_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t):
    
    p=p_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t)
    alpha=3
    rho=rho_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t)
    u=u_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t)
    
    return 0.5*alpha*p + 0.5*rho*u**2


def plot_E_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t):
    
    plt.plot(x,E_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t),lw=3)

    Emax=max(E_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t))
    Emin=min(E_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t))
    Ediff=Emax-Emin
    fts=30
    plt.axis([xmin,xmax,Emin-0.1*Ediff,Emax+0.1*Ediff])
    plt.xlabel(r'$x$',fontsize=fts)
    plt.ylabel(r'$E$',fontsize=fts)
    plt.title(r'$t=%.3f$' % t,fontsize=fts)

def plot_rho_u_p_E_iso(rhol,ul,pl,cl,rhor,ur,pr,t):
    
    plt.figure(figsize=(26,5))

    plt.subplot(141)
    plot_rho_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t)

    plt.subplot(142)
    plot_u_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t)
    
    plt.subplot(143)
    plot_p_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t)

    plt.subplot(144)
    plot_E_lin_iso(rhol,ul,pl,cl,rhor,ur,pr,t)

plot_rho_u_p_E_iso(1,0,1,1,2,0,2,0.01)




