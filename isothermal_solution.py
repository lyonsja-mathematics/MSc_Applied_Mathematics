# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 10:54:38 2016

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

def rie_1d_isoth_exact(rhol,ul,rhor,ur,T,t,x):
    """Exact solution to 1d Isentropic equations
    """
    
    
    #obtain isentropic constant for internal use
    #kap=cl**2/(gamma*rhol**(gamma-1))
    R=8.314
    
    def speed(T):
        return np.sqrt(R*T)
    def press(rho,T):
        return rho*R*T
    
    pl=press(rhol,T)
    pr=press(rhor,T)
    a=speed(T)
    
    
    #expressions for downstream velocity ustar as function of downstream density rhostar
    shk1 = lambda rhostar : ul - np.sqrt((1/rhol-1/rhostar)*(press(rhostar,T)-pl))
    rar1 = lambda rhostar : ul -a*np.log(rhostar/rhol) 
    
    shk2 = lambda rhostar : ur + np.sqrt((1/rhor-1/rhostar)*(press(rhostar,T)-pr))
    rar2 = lambda rhostar : ur +a*np.log(rhostar/rhor)
     
    
    #decide on whether wave is shock or rarefaction for given value of rhostar
    def phil(rhostar):        
        if rhostar<rhol: return rar1(rhostar)
        else: return shk1(rhostar)
    
    def phir(rhostar):
        if rhostar<rhor: return rar2(rhostar)
        else: return shk2(rhostar)
        
    #difference in values of downstream velocities is measure of error in rhostar    
    phi = lambda rhostar : phil(rhostar)-phir(rhostar)

    #rhostar=1.4
    rhostar,info, ier, msg = fsolve(phi, (rhol+rhor)/2.,full_output=True,factor=0.1,xtol=1.e-10)
    
    #common upstream velocity
    ustar = phil(rhostar)
    
    #common upstream soundspeed
    cstar=speed(rhostar)
    
    ws = np.zeros(4) 
    
    if rhostar<rhol: 
        ws[0] = ul - a
        ws[1] = ustar - a
    else:
        ws[0] = (rhol*ul - rhostar*ustar)/(rhol - rhostar)
        ws[1] = ws[0]
  
    if rhostar<rhor: 
        ws[2] = ustar + a
        ws[3] = ur + a
    else:
        ws[2] = (rhor*ur - rhostar*ustar)/(rhor - rhostar)
        ws[3] = ws[2]
  

    xs = ws*t # Wave front positions
    
    xi = x/t  # coordinate for rarefactions
    
    u1=a+xi                                             #left moving rarefaction
    #c1=cl+0.5*(gamma-1)*(ul-u1)
    rho1=rhol*np.exp((ul-u1)/a)

    u2=xi-a                                             #right moving rarefaction
    #c2=cr-0.5*(gamma-1)*(ur-u2)
    rho2=rhor*np.exp((u2-ur)/a)
    
    cnd=[x<=xs[0],(x>xs[0]) & (x<xs[1]),(x>=xs[1]) & (x<=xs[2]),(x>xs[2]) & (x<xs[3]),x>=xs[3]]
    
    uout=np.select(cnd,[ul,u1,ustar,u2,ur])
    #cout=np.select(cnd,[cl,c1,cstar,c2,cr])
    rhoout=np.select(cnd,[rhol,rho1,rhostar,rho2,rhor])

    
    return rhoout, uout, rhostar
    
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
    
def plot_rie_1d_euler_exact(rhol,ul,rhor,ur,T,t):

    x = np.linspace(-1.,1.,1000)
    
    rho, u, rhostar = rie_1d_isoth_exact(rhol,ul,rhor,ur,T,t,x)

    rholin= rho_lin(rhol,ul,rhor,ur,t)
    ulin=u_lin(rhol,ul,rhor,ur,t)
    
    fig = plt.figure(figsize=(15,6))
    sol = [rho,u]
    linsol=[rholin,ulin]
    names = ['Density','Velocity']
    ylabels = [r'$\rho$','$u$']
    axes = [0]*3
    fts=30
    for i in range(2):
        axes[i] = fig.add_subplot(1,2,i+1)
        q = sol[i]
        qlin=linsol[i]
        plt.plot(x,q,x,qlin,linewidth=3)
        plt.title(r'$t=%.3f$' % t, fontsize=fts)
        qmax = max(q)
        qmin = min(q)
        qdiff = qmax - qmin
        axes[i].set_ylim((qmin-0.1*qdiff,qmax+0.1*qdiff))
        plt.xlabel(r'$x$',fontsize=fts)
        plt.ylabel(ylabels[i],fontsize=fts)
        plt.legend((r'Non Linear','Linear'),loc='best',fontsize=0.5*fts)
    #print(rhostar)






#for t in [0,0.01,0.05]:
plot_rie_1d_euler_exact(1,0,2,0,293,0.02)