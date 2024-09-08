
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

def rie_1d_isen_exact(rhol,ul,pl,rhor,ur,pr,gamma,t,x):
    """Exact solution to 1d Isentropic equations
    """
    
    
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
        return  p/(gamma-1.) + 0.5*rho*u**2

    El=E(rhol,ul,pl,gamma)
    Elstar=E(rholstar,ustar,pstar,gamma)
    Erstar=E(rhorstar,ustar,pstar,gamma)
    Er=E(rhor,ur,pr,gamma)
    
    u1=((gamma-1)/(gamma+1))*ul+(2/(gamma+1))*(xi+cl) #left moving rarefaction
    c1=cl+0.5*(gamma-1)*(ul-u1)
    rho1=rhol*(c1/cl)**(2/(gamma-1))
    p1=(1/gamma)*rho1*c1**2
    E1=E(rho1,u1,p1,gamma)

    u2=((gamma-1)/(gamma+1))*ur+(2/(gamma+1))*(xi-cr) #right moving rarefaction
    c2=cr-0.5*(gamma-1)*(ur-u2)
    rho2=rhor*(c2/cr)**(2/(gamma-1))
    p2=(1/gamma)*rho2*c2**2
    E2=E(rho2,u2,p2,gamma)
    
    cnd=[x<=xs[0], (x>xs[0]) & (x<xs[1]), (x>=xs[1]) & (x<=xs[2]), (x>xs[2]) & (x<xs[3]), (x>xs[3]) & (x<xs[4]), x>=xs[4]]
    
    
    

    uout=np.select(cnd,[ul,u1,ustar,ustar,u2,ur])
    cout=np.select(cnd,[cl,c1,clstar,crstar,c2,cr])
    rhoout=np.select(cnd,[rhol,rho1,rholstar,rhorstar,rho2,rhor])
    pout=np.select(cnd,[pl,p1,pstar,pstar,p2,pr])
    #Eout=np.select(cnd,[El,E1,Elstar,Erstar,E2,Er])
    Eout=E(rhoout,uout,pout,gamma)
    
    

    
    return rhoout, uout, cout, pout, Eout, pstar, ustar
    
def plot_rie_1d_euler_exact(rhol,ul,pl,rhor,ur,pr,gamma,t):

    x = np.linspace(-1.,1.,1000)
    
    rho, u, c, p, E, pstar, ustar = rie_1d_isen_exact(rhol,ul,pl,rhor,ur,pr,gamma,t,x)

    fig = plt.figure(figsize=(7,5))
    sol = [rho,u,p,E]
    #names = ['Density','Velocity','Pressure','Energy']
    ylabels = [r'$\rho$','$u$','$q$','$E$']
    axes = [0]*4
    fts=30
    #for i in range(2,3):
        #axes[i]=fig
        #axes[i] = fig.add_subplot(1,1,i+1)
    q = sol[2]
    plt.plot(x,q,linewidth=3)
        #plt.title(r'$t=%.2f$' % t, fontsize=fts)
    qmax = max(q)
    qmin = min(q)
    qdiff = qmax - qmin
    plt.axis([-1,1,qmin-0.1*qdiff,qmax+0.1*qdiff])
        #axes[i].set_ylim((qmin-0.1*qdiff,qmax+0.1*qdiff))
    plt.xlabel(r'$x$',fontsize=fts)
    plt.ylabel(ylabels[2],fontsize=fts)
    
    plt.suptitle(r'$t>0$', fontsize=fts)
    print(pstar,ustar)
        

#rar_rar=plot_rie_1d_euler_exact(1,-1,1,1,1,1,1.4,0.4)
#shk_rar=plot_rie_1d_euler_exact(1,0,1,3,0,3,1.4,0.4)
rar_shk=plot_rie_1d_euler_exact(3,0,3,1,0,1,1.4,0.4)
#shk_shk=plot_rie_1d_euler_exact(1,1,1,1,-1,1,1.4,0.4)














