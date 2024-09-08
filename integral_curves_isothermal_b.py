# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 14:33:35 2016

@author: lyons
"""

get_ipython().magic('matplotlib inline')
from numpy import linspace,zeros,exp,arange,sqrt, cos,log
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
def plot_Hugoniot_locus(rhol,ul,rhor,ur):
    
    R=8.314
    T=293
    a=sqrt(R*T)
    
    
    
    def p(rho):
        return (a**2)*rho
    
    
    #expressions for downstream velocity ustar as function of downstream density rhostar
    shk1 = lambda rhostar : ul - sqrt((((1/rhol)-(1/rhostar))*(p(rhostar)-p(rhol))))
    rar1 = lambda rhostar : ul - a*log(rhostar/rhol)
    
    shk3 = lambda rhostar : ur + sqrt((((1/rhor)-(1/rhostar))*(p(rhostar)-p(rhor))))
    rar3 = lambda rhostar : ur + a*log(rhostar/rhor)
     
    
    #decide on whether wave is shock or rarefaction for given value of rhostar
    def phil(rhostar):        
        if rhostar<rhol: return rar1(rhostar)
        else: return shk1(rhostar)
    
    def phir(rhostar):
        if rhostar<rhor: return rar3(rhostar)
        else: return shk3(rhostar)
        
    #difference in values of downstream velocities is measure of error in rhostar    
    phi = lambda rhostar : phil(rhostar)-phir(rhostar)

    rhostar,info, ier, msg = fsolve(phi, (rhol+rhor)/2.,full_output=True,factor=0.1,xtol=1.e-10)
    
    #common upstream velocity
    ustar = phil(rhostar)
    
    
    figure(figsize=(7,6))
    
    rho=np.linspace(0,4,1000)
    
    plot(rho,shk1(rho),rho,rar3(rho),linewidth=2)
    plot(rhostar,ustar,marker='o',ms=12,color='black')
    plot(rhol,ul,marker='o',ms=12,color='blue')
    plot(rhor,ur,marker='o',ms=12,color='green')
    axis([0,3,-100,50])
    fts=30
    xlabel(r'$\rho^*$', fontsize=fts)
    ylabel(r'$u^*$', fontsize=fts)
    print(rhostar,ustar)
    legend((r'1-shock locus','2-rarefaction locus'),fontsize=fts*0.5,loc='best')
    annotate(r'$(1.41,-17.15)$', xy=(1.1,-40),fontsize=fts*0.6)
    annotate(r'$(1,0)$', xy=(0.9,10),fontsize=fts*0.6)
    annotate(r'$(2,0)$', xy=(2,-15),fontsize=fts*0.6)
    #annotate(r'$(1,0)$', xy=(1,-0.3),fontsize=20)
    #annotate(r'$\left( \frac{1}{2} , 0\right)$', xy=(0.5,-0.3),fontsize=20)
    
    
plot_Hugoniot_locus(1,0,2,0)