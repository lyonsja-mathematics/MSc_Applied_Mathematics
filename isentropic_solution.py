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

def rie_1d_isen_exact(rhol,ul,cl,rhor,ur,gamma,t,x):
    """Exact solution to 1d Isentropic equations
    """
    
    
    #obtain isentropic constant for internal use
    kap=cl**2/(gamma*rhol**(gamma-1))
    
    def speed(rho):
        return np.sqrt(kap*gamma*rho**(gamma-1))
    def press(rho):
        return kap*rho**gamma
    
    pl=press(rhol)
    pr=press(rhor)
    cr=speed(rhor)
    
    
    #expressions for downstream velocity ustar as function of downstream density rhostar
    shk1 = lambda rhostar : ul - np.sqrt((1/rhostar-1/rhol)*(pl-press(rhostar)))
    rar1 = lambda rhostar : ul +(2/(gamma-1))*(cl-speed(rhostar)) 
    
    shk2 = lambda rhostar : ur + np.sqrt((1/rhostar-1/rhor)*(pr-press(rhostar)))
    rar2 = lambda rhostar : ur -(2/(gamma-1))*(cr-speed(rhostar)) 
     
    
    #decide on whether wave is shock or rarefaction for given value of rhostar
    def phil(rhostar):        
        if rhostar<rhol: return rar1(rhostar)
        else: return shk1(rhostar)
    
    def phir(rhostar):
        if rhostar<rhor: return rar2(rhostar)
        else: return shk2(rhostar)
        
    #difference in values of downstream velocities is measure of error in rhostar    
    phi = lambda rhostar : phil(rhostar)-phir(rhostar)

    rhostar,info, ier, msg = fsolve(phi, (rhol+rhor)/2.,full_output=True,factor=0.1,xtol=1.e-10)
    
    #common upstream velocity
    ustar = phil(rhostar)
    
    #common upstream soundspeed
    cstar=speed(rhostar)
    
    ws = np.zeros(4) 
    
    if rhostar<rhol: 
        ws[0] = ul - cl
        ws[1] = ustar - cstar
    else:
        ws[0] = (rhol*ul - rhostar*ustar)/(rhol - rhostar)
        ws[1] = ws[0]
  
    if rhostar<rhor: 
        ws[2] = ustar + cstar
        ws[3] = ur + cr
    else:
        ws[2] = (rhor*ur - rhostar*ustar)/(rhor - rhostar)
        ws[3] = ws[2]
  

    xs = ws*t # Wave front positions
    
    xi = x/t  # coordinate for rarefactions
    
    u1=((gamma-1)/(gamma+1))*ul+(2/(gamma+1))*(xi+cl) #left moving rarefaction
    c1=cl+0.5*(gamma-1)*(ul-u1)
    rho1=rhol*(c1/cl)**(2/(gamma-1))

    u2=((gamma-1)/(gamma+1))*ur+(2/(gamma+1))*(xi-cr) #right moving rarefaction
    c2=cr-0.5*(gamma-1)*(ur-u2)
    rho2=rhol*(c2/cl)**(2/(gamma-1))
    
    cnd=[x<=xs[0],(x>xs[0]) & (x<xs[1]),(x>=xs[1]) & (x<=xs[2]),(x>xs[2]) & (x<xs[3]),x>=xs[3]]
    
    uout=np.select(cnd,[ul,u1,ustar,u2,ur])
    cout=np.select(cnd,[cl,c1,cstar,c2,cr])
    rhoout=np.select(cnd,[rhol,rho1,rhostar,rho2,rhor])
    pout=press(rhoout)

    
    return rhoout, uout, cout, pout
    
def plot_rie_1d_euler_exact(rhol,ul,cl,rhor,ur,gamma,t):

    x = np.linspace(-1.,1.,1000)
    
    rho, u, c, p = rie_1d_isen_exact(rhol,ul,cl,rhor,ur,gamma,t,x)

    fig = plt.figure(figsize=(20,6))
    sol = [rho,u,c]
    names = ['Density','Velocity','Soundspeed']
    ylabels = [r'\rho','u','c']
    axes = [0]*3
    for i in range(3):
        axes[i] = fig.add_subplot(1,3,i+1)
        q = sol[i]
        plt.plot(x,q,linewidth=3)
        plt.title(r'$t=%.2f$' % t, fontsize=20)
        qmax = max(q)
        qmin = min(q)
        qdiff = qmax - qmin
        axes[i].set_ylim((qmin-0.1*qdiff,qmax+0.1*qdiff))
        plt.xlabel(r'x',fontsize=20)
        plt.ylabel(ylabels[i],fontsize=20)

    #print(rhostar,ustar,cstar)



#for t in np.arange(0,1.25,0.25):
#plot_rie_1d_euler_exact(1,0,1,0.1,0,1.4,0.33)
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

'''
plt.subplot(223)
plt.plot(x,c,lw=3)
plt.xlabel(r'$x$',fontsize=25)
plt.ylabel(r'$c$',fontsize=25)
plt.axis([-1,1,0.6,1.1])
'''

plt.subplot(133)
plt.plot(x,p,lw=3)
plt.xlabel(r'$x$',fontsize=fts)
plt.ylabel(r'$p$',fontsize=fts)









