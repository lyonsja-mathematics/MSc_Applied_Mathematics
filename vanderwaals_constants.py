# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 11:57:02 2016

@author: lyons
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def plot_ab():
    
    bstart=0
    bend=0.2
    npts=1000
    b=np.linspace(bstart,bend,npts)
    
    def a(b,T,rho):
        return (b*8.314*T)/(1-b*rho)
    
    T_list=[273,293,393]
    rho_list=[0.5,1,2]
      
    fig=plt.figure(figsize=(10,4))
    axes=[0]*2

    fts=20
    
    axes[0]=fig.add_subplot(121)
    for T in T_list:
            plt.plot(b,a(b,T,1),lw=2)
            ymax=3
            ymin=0
            ydiff=ymax-ymin
            axes[0].set_ylim((ymin,ymax))
            axes[0].set_xlim((bstart,bend))
            plt.xlabel(r'b',fontsize=fts)
            plt.ylabel(r'a',fontsize=fts)
            plt.legend((r'T=273 K','T=293 K','T=393 K'),loc='best',fontsize=0.5*fts)
            
        
    axes[1]=fig.add_subplot(122)
    for rho in rho_list:
        plt.plot(b,a(b,293,rho))
        
plot_ab()