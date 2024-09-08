
# coding: utf-8

# In[14]:

get_ipython().magic('matplotlib inline')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math as m
import warnings
warnings.filterwarnings("ignore")

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def rie_1d_isen_lin(rhol,ul,pl,rhor,ur,pr,gamma,t,x):
    """Exact solution to 1d Isentropic equations
    """
    
    npt=1000
    #obtain isentropic constant for internal use
    kapl=pl/(rhol**gamma)
    kapr=pr/(rhor**gamma)
    
    def speed_l(p):
        return np.sqrt(gamma*kapl**(1/gamma))*p**((gamma-1)/(2*gamma))
    def press_l(rho):
        return kapl*rho**gamma
    def rho_l(p):
        return (p/kapl)**(1/gamma)
        
    def speed_r(p):
        return np.sqrt(gamma*kapr**(1/gamma))*p**((gamma-1)/(2*gamma))
    def press_r(rho):
        return kapr*rho**gamma
    def rho_r(p):
        return (p/kapr)**(1/gamma)
    
    #pl=press_l(rhol)
    cl=speed_l(pl)
    #pr=press_r(rhor)
    cr=speed_r(pr)
    
    
    #expressions for downstream velocity ustar as function of downstream density rhostar
    shk1 = lambda pstar : ul - np.sqrt((((1/rho_l(pstar))-(1/rhol))*(pl-pstar)))
    rar1 = lambda pstar : ul + cl - (cl/rhol)*rho_l(pstar)
    
    shk3 = lambda pstar : ur + np.sqrt((((1/rho_r(pstar))-(1/rhor))*(pr-pstar)))
    rar3 = lambda pstar : ur - cr + (cr/rhor)*rho_r(pstar) 
     
    
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
    
    def H(x):
        v=np.zeros(npt)
        for i in range(npt):
            if x[i]<0:
                v[i]=0
            if x[i]>0:
                v[i]=1
        return v
        
    def w_l(rhol,ul,pl,cl,rhor,ur,pr):
        u0=ustar
        rho0=(rhorstar+rholstar)/2
        c0=(crstar+clstar)/2
        
        return rhol-(pl/(c0**2)), 0.5*ul - (pl/(2*rho0*c0)), 0.5*ul + (pl/(2*rho0*c0))
        
                
                
    def w_r(rhol,ul,pl,cl,rhor,ur,pr):
        u0=ustar
        rho0=(rhorstar+rholstar)/2
        c0=(crstar+clstar)/2
        
        return rhor-(pr/(c0**2)), 0.5*ur - (pr/(2*rho0*c0)), 0.5*ur + (pr/(2*rho0*c0))
       
    def u_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t):
        
        u0=ustar
        rho0=(rhorstar+rholstar)/2
        c0=(crstar+clstar)/2
        
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

  
    def rho_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t):
        
        u0=ustar
        rho0=(rhorstar+rholstar)/2
        c0=(crstar+clstar)/2
        
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
        
    
        
    def p_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t):
        
        u0=ustar
        rho0=(rhorstar+rholstar)/2
        c0=(crstar+clstar)/2
        
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
    
   
    def E_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t):
        
        p=p_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t)
        gamma=1.4
        rho=rho_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t)
        u=u_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t)
        
        return (p/(gamma-1)) + 0.5*rho*u**2

    rhoout=rho_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t)
    uout=u_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t)
    pout=p_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t)
    Eout=E_lin_isen(rhol,ul,pl,cl,rhor,ur,pr,t)
    
    return rhoout,uout,pout,Eout
    
    

    
    return rhoout, uout, pout, Eout
    
def plot_rie_1d_euler_lin(rhol,ul,pl,rhor,ur,pr,gamma,t):

    x = np.linspace(-1.,1.,1000)
    
    rho, u, p, E = rie_1d_isen_lin(rhol,ul,pl,rhor,ur,pr,gamma,t,x)

    fig = plt.figure(figsize=(20,5))
    sol = [rho,u,p,E]
    #names = ['Density','Velocity','Pressure','Energy']
    ylabels = [r'$\rho$','$u$','$p$','$E$']
    axes = [0]*4
    fts=30
    for i in range(4):
        axes[i] = fig.add_subplot(1,4,i+1)
        q = sol[i]
        plt.plot(x,q,linewidth=3)
        #plt.title(r'$t=%.2f$' % t, fontsize=fts)
        qmax = max(q)
        qmin = min(q)
        qdiff = qmax - qmin
        axes[i].set_ylim((qmin-0.1*qdiff,qmax+0.1*qdiff))
        plt.xlabel(r'$x$',fontsize=fts)
        plt.ylabel(ylabels[i],fontsize=fts)
    
    plt.suptitle(r'$t=%.2f$' % t, fontsize=fts)
    
        

plot_rie_1d_euler_lin(3,0,3,1,0,1,1.4,0.4)