# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 16:56:20 2016

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

def rie_1d_isen_exact(rhol,ul,cl,rhor,ur,cr,gamma,t,x):
    """Exact solution to 1d Isentropic equations
    """
    
    
    #obtain isentropic constant for internal use
    kapl=cl**2/(gamma*rhol**(gamma-1))
    kapr=cr**2/(gamma*rhor**(gamma-1))
    
    def speed_l(rho):
        return np.sqrt(kapl*gamma*rho**(gamma-1))
    def press_l(rho):
        kapl*rho**gamma
    def rho_l(p):
        return (p/kapl)**(1/gamma)
        
    def speed_r(rho):
        return np.sqrt(kapr*gamma*rho**(gamma-1))
    def press_r(rho):
        kapr*rho**gamma
    def rho_r(p):
        return (p/kapr)**(1/gamma)
    
    pl=press_l(rhol)
    #cl=speed_l(rhol)
    pr=press_r(rhor)
    #cr=speed_r(rhor)
    
    
    #expressions for downstream velocity ustar as function of downstream density rhostar
    shk1 = lambda pstar : ul - np.sqrt((((1/rho_l(pstar))-(1/rhol))*(pl-pstar)))
    rar1 = lambda pstar : ul + (2*cl/(gamma-1))*(1-(pstar/pl)**((gamma-1)/(2*gamma)))
    
    shk3 = lambda pstar : ur + np.sqrt((((1/rho_r(pstar))-(1/rhor))*(pl-pstar)))
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
    
    
    #def f(rholstar):
     #   return press_l(rholstar)-pstar
    rholstar=rho_l(pstar)
    
    
    #g = lambda rhorstar : press_r(rhorstar)-pstar
    rhorstar = rho_r(pstar)
    
    #common upstream soundspeed
    clstar=speed_l(pstar)
    crstar=speed_r(pstar)
    
    ws = np.zeros(5) 
    
    if pstar<pl: 
        ws[0] = ul - cl
        ws[1] = ustar - clstar
    else:
        ws[0] = (rhol*ul - rholstar*ustar)/(rhol - rholstar)
        ws[1] = ws[0]

    ws[2]=ustar
  
    if pstar<pr: 
        ws[3] = ustar + crstar
        ws[4] = ur + cr
    else:
        ws[3] = (rhor*ur - rhorstar*ustar)/(rhor - rhorstar)
        ws[4] = ws[3]
  

    xs = ws*t # Wave front positions
    
    xi = x/t  # coordinate for rarefactions
    
    def E(rho,u,p,gamma):
        return  0.5*rho*u**2

    El=E(rhol,ul,pl,gamma)
    Elstar=E(rholstar,ustar,pstar,gamma)
    Erstar=E(rhorstar,ustar,pstar,gamma)
    Er=E(rhor,ur,pr,gamma)
    
    u1=((gamma-1)/(gamma+1))*ul+(2/(gamma+1))*(xi+cl) #left moving rarefaction
    c1=cl+0.5*(gamma-1)*(ul-u1)
    rho1=rhol*(c1/cl)**(2/(gamma-1))
    p1=press_l(rho1)
    E1=E(rho1,u1,p1,gamma)

    u2=((gamma-1)/(gamma+1))*ur+(2/(gamma+1))*(xi-cr) #right moving rarefaction
    c2=cr-0.5*(gamma-1)*(ur-u2)
    rho2=rhol*(c2/cr)**(2/(gamma-1))
    p2=press_r(rho2)
    E2=E(rho2,u2,p2,gamma)
    
    cnd=[x<=xs[0], (x>xs[0]) & (x<xs[1]), (x>=xs[1]) & (x<=xs[2]), (x>xs[2]) & (x<=xs[3]), (x>xs[3]) & (x<xs[4]), x>=xs[4]]
    
    
    

    uout=np.select(cnd,[ul,u1,ustar,ustar,u2,ur])
    cout=np.select(cnd,[cl,c1,clstar,crstar,c2,cr])
    rhoout=np.select(cnd,[rhol,rho1,rholstar,rhorstar,rho2,rhor])
    pout=np.select(cnd,[pl,p1,pstar,pstar,p2,pr])
    Eout=np.select(cnd,[El,E1,Elstar,Erstar,E2,Er])
    
    

    
    return rhoout, uout, cout, pout, Eout
    
def plot_rie_1d_euler_exact(rhol,ul,cl,rhor,ur,cr,gamma,t):

    x = np.linspace(-1.,1.,1000)
    
    rho, u, c, p, E = rie_1d_isen_exact(rhol,ul,cl,rhor,ur,cr,gamma,t,x)

    fig = plt.figure(figsize=(15,12))
    sol = [rho,u,p,E]
    #names = ['Density','Velocity','Pressure','Energy']
    ylabels = [r'$\rho$','$u$','$p$','$E$']
    axes = [0]*4
    for i in range(4):
        axes[i] = fig.add_subplot(2,2,i+1)
        q = sol[i]
        plt.plot(x,q,linewidth=3)
        plt.title(r'$t=%.2f$' % t, fontsize=20)
        qmax = max(q)
        qmin = min(q)
        qdiff = qmax - qmin
        axes[i].set_ylim((qmin-0.1*qdiff,qmax+0.1*qdiff))
        plt.xlabel(r'$x$',fontsize=20)
        plt.ylabel(ylabels[i],fontsize=20)


plot_rie_1d_euler_exact(1,0,1,0.1,0,0.5,1.4,0)
'''
x=x = np.linspace(-1.,1.,1000)
rho, u, c, p = rie_1d_isen_exact(1,0,1,0.1,0,1.4,0.4,x)

plt.figure(figsize=(20,6))

fts=30
#plt.title(r'$t=0.33$',fontsize=25)

plt.subplot(131)
plt.plot(x,rho,lw=3)
plt.xlabel(r'$x$',fontsize=fts)
plt.ylabel(r'$\rho$',fontsize=fts)
plt.axis([-1,1,0,1.1])

plt.subplot(132)
plt.plot(x,u,lw=3)
plt.xlabel(r'$x$',fontsize=fts)
plt.ylabel(r'$u$',fontsize=fts)


plt.subplot(223)
plt.plot(x,c,lw=3)
plt.xlabel(r'$x$',fontsize=25)
plt.ylabel(r'$c$',fontsize=25)
plt.axis([-1,1,0.6,1.1])

plt.subplot(133)
plt.plot(x,p,lw=3)
plt.xlabel(r'$x$',fontsize=fts)
plt.ylabel(r'$p$',fontsize=fts)
'''