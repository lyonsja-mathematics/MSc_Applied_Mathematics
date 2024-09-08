# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 14:33:35 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
from numpy import linspace,zeros,exp,arange,sqrt, cos
from warnings import filterwarnings
from math import pi
from scipy.optimize import fsolve
from matplotlib.pyplot import *
filterwarnings("ignore")
rc('text',usetex=True)
rc('font', family='serif')


'''   
def v_2(rho_0,u_0,gamma,c_0,rho):
    return u_0 + (2*c_0/(gamma-1))*(1-(rho/rho_0)**((gamma-1)/(2*gamma)))
    
def v_3(rho_0,u_0,gamma,c_0,rho):
    return u_0 - (2*c_0/(gamma-1))*(1-(rho/rho_0)**((gamma-1)/(2*gamma)))
'''
def plot_Hugoniot_locus(rhol,ul,pl,rhor,ur,pr):
    
    gamma=1.4
    kapl=pl/(rhol**gamma)
    kapr=pr/(rhor**gamma)
    
    
    
    def press_l(rho):
        return kapl*rho**gamma
    def rho_l(p):
        return (p/kapl)**(1/gamma)
        
    
    def press_r(rho):
        return kapr*rho**gamma
    def rho_r(p):
        return (p/kapr)**(1/gamma)
    
    #pl=press_l(rhol)
    cl=np.sqrt(gamma*pl/rhol)
    #pr=press_r(rhor)
    cr=np.sqrt(gamma*pl/rhol)
    
    
    #expressions for downstream velocity ustar as function of downstream density rhostar
    shk1 = lambda pstar : ul - np.sqrt((((1/rho_l(pstar))-(1/rhol))*(pl-pstar)))
    rar1 = lambda pstar : ul + (2*cl/(gamma-1))*(1-(pstar/pl)**((gamma-1)/(2*gamma)))
    
    shk3 = lambda pstar : ur + np.sqrt((((1/rho_r(pstar))-(1/rhor))*(pr-pstar)))
    rar3 = lambda pstar : ur - (2*cr/(gamma-1))*(1-(pstar/pr)**((gamma-1)/(2*gamma))) 
     
    
    #decide on whether wave is shock or rarefaction for given value of rhostar
    def phil(pstar):        
        if pstar<pl: return rar1(pstar)
        else: return shk1(pstar)
    
    def phir(pstar):
        if pstar<pr: return rar3(pstar)
        else: return shk3(pstar)
        
    #difference in values of downstream velocities is measure of error in rhostar    
    phi = lambda pstar : phil(pstar)-phir(pstar)

    pstar,info, ier, msg = fsolve(phi, (pl+pr)/2.,full_output=True,factor=0.1,xtol=1.e-10)
    
    #common upstream velocity
    ustar = phil(pstar)
    
    
    plt.figure(figsize=(10,6))
    
    p=np.linspace(0,4,1000)
    
    plt.plot(p,rar1(p),p,shk3(p),linewidth=3)
    plt.plot(pstar,ustar,marker='o',ms=12,color='black')
    plt.plot(pl,ul,marker='o',ms=12,color='blue')
    plt.plot(pr,ur,marker='o',ms=12,color='green')
    axis([0,4,-0.25,1])
    fts=30
    plt.xlabel(r'$p$', fontsize=fts)
    plt.ylabel(r'$u$', fontsize=fts)
    #print(pstar,ustar)
    legend((r'2-rarefaction locus','3-shock locus'),fontsize=fts*0.5,loc='best')
    annotate(r'$(1.69,0.47)$', xy=(1.9,0.45),fontsize=fts*0.6)
    annotate(r'$(3,0)$', xy=(3,0.1),fontsize=fts*0.6)
    annotate(r'$(1,0)$', xy=(1.1,0),fontsize=fts*0.6)
    #annotate(r'$(1,0)$', xy=(1,-0.3),fontsize=20)
    #annotate(r'$\left( \frac{1}{2} , 0\right)$', xy=(0.5,-0.3),fontsize=20)
    
    
plot_Hugoniot_locus(3,0,3,1,0,1)