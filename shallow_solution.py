# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 11:49:42 2016

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

def rie_1d_shall_exact(hl,ul,hr,ur,t,x):
    """Exact solution to 1d Isentropic equations
    """
    g=9.81
    rho=1
    
    #obtain isentropic constant for internal use
    #kap=cl**2/(gamma*rhol**(gamma-1))
    
    def speed(h):
        return np.sqrt(g*h)
    def press(h):
        return rho*g*h
    
    pl=press(hl)
    pr=press(hr)
    
    cl=speed(hl)
    cr=speed(hr)
    
    def u_star(hl,ul,hr,ur):
        return (ul+ur)/2 +np.sqrt(hl/g)-np.sqrt(hr/g)
        
    def h_star(hl,ul,hr,ur):
        return (np.sqrt(hr)+np.sqrt(hl))/2 +(ul-ur)/(4*np.sqrt(g))
    
    #ustar=u_star(hl,ul,hr,ur)
    #hstar=h_star(hl,ul,hr,ur)
    #cstar=speed(ustar)
    #expressions for downstream velocity ustar as function of downstream density rhostar
    shk1 = lambda hstar : ul-np.sqrt(g/2*((hstar/hl)-(hl/hstar))*(hstar-hl)) 
    rar1 = lambda hstar : ul - 2*np.sqrt(g)*(np.sqrt(hstar)-np.sqrt(hl))
    
    shk2 = lambda hstar : ur+np.sqrt(g/2*((hstar/hr)-(hr/hstar))*(hstar-hr))
    rar2 = lambda hstar : ur + 2*np.sqrt(g)*(np.sqrt(hstar)-np.sqrt(hr))
     
    
    #decide on whether wave is shock or rarefaction for given value of rhostar
    def phil(hstar):        
        if hstar<hl: return rar1(hstar)
        else: return shk1(hstar)
    
    def phir(hstar):
        if hstar<hr: return rar2(hstar)
        else: return shk2(hstar)
        
    #difference in values of downstream velocities is measure of error in rhostar    
    phi = lambda hstar : phil(hstar)-phir(hstar)

    hstar,info, ier, msg = fsolve(phi, (hl+hr)/2.,full_output=True,factor=0.1,xtol=1.e-10)
    
    #common upstream velocity
    ustar = phil(hstar)
    
    #common upstream soundspeed
    cstar=speed(hstar)
    
    ws = np.zeros(4) 
    
    if hstar<hl:
        ws[0] = ul - cl
        ws[1] = ustar - cstar
    else:
        ws[0] = (hl*ul - hstar*ustar)/(hl - hstar)
        ws[1] = ws[0]
  
    if hstar<hr: 
        ws[2] = ustar + cstar
        ws[3] = ur + cr
    else:
        ws[2] = (hr*ur - hstar*ustar)/(hr - hstar)
        ws[3] = ws[2]
  

    xs = ws*t # Wave front positions
    
    xi = x/t  # coordinate for rarefactions
    
    #left moving rarefaction
    u1=(ul/3)+(2/3)*(xi+np.sqrt(g*hl)) 
    #c1=cl+0.5*(gamma-1)*(ul-u1)
    h1=(1/g)*(u1-xi)**2
    
    #right moving rarefaction
    u2=(ur/3)+(2/3)*(xi-np.sqrt(g*hr))  
    #c2=cr-0.5*(gamma-1)*(ur-u2)
    h2=(1/g)*(u2-xi)**2
    
    cnd=[x<=xs[0],(x>xs[0]) & (x<xs[1]),(x>=xs[1]) & (x<=xs[2]),(x>xs[2]) & (x<xs[3]),x>=xs[3]]
    
    uout=np.select(cnd,[ul,u1,ustar,u2,ur])
    #cout=np.select(cnd,[cl,c1,cstar,c2,cr])
    hout=np.select(cnd,[hl,h1,hstar,h2,hr])
    pout=press(hout)

    
    return hout, uout, pout
    
def plot_rie_1d_euler_exact(hl,ul,hr,ur,t):

    x = np.linspace(-2.,2.,2000)
    '''
    def h(t):
        return rie_1d_shall_exact(hl,ul,hr,ur,t,x)[0]

    def u(t):
        return rie_1d_shall_exact(hl,ul,hr,ur,t,x)[1]

    def p(t):
        return rie_1d_shall_exact(hl,ul,hr,ur,t,x)[2]
        
    #fig = plt.figure(figsize=(20,6))
    '''
    h, u, p = rie_1d_shall_exact(hl,ul,hr,ur,t,x)

    fig = plt.figure(figsize=(20,6))
    sol = [h,u,p]
    #names = ['Density','Velocity','Soundspeed']
    ylabels = [r'$h$','$u$','$p$']
    axes = [0]*3
    for i in range(3):
        axes[i] = fig.add_subplot(1,3,i+1)
        q = sol[i]
        plt.plot(x,q,linewidth=3)
        plt.title(r'$t=%.2f$' % t, fontsize=25)
        qmax = max(q)
        qmin = min(q)
        qdiff = qmax - qmin
        axes[i].set_ylim((qmin-0.1*qdiff,qmax+0.1*qdiff))
        plt.xlabel(r'$x$',fontsize=25)
        plt.ylabel(ylabels[i],fontsize=25)

    #print(rhostar,ustar,cstar)



#for t in np.arange(0,1.25,0.25):
plot_rie_1d_euler_exact(0.5,-1,0.5,1,0.5)
'''
def plot_b(hl,ul,hr,ur,t):
    
    x = np.linspace(-2.,2.,2000)
    
    def h(t):
        return rie_1d_shall_exact(hl,ul,hr,ur,t,x)[0]

    def u(t):
        return rie_1d_shall_exact(hl,ul,hr,ur,t,x)[1]

    def p(t):
        return rie_1d_shall_exact(hl,ul,hr,ur,t,x)[2]
        
    fig = plt.figure(figsize=(20,6))
    
    plt.subplot(131)
    for t in [0,0.3,0.5,1]:
        plt.plot(x, h(t),lw=3)
    plt.legend((r'$t=%.2f$' % t), fontsize=12,loc='best')
    

plot_b(0.5,-1,0.5,1,0)
'''   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    