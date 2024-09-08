# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 19:11:09 2016

@author: lyons
"""

from matplotlib.pyplot import *
from numpy import linspace, sqrt,cos,sin,log
import warnings
warnings.filterwarnings("ignore")

rc('text', usetex=True)
rc('font', family='serif')

dx=0.02*pi; r=0.06
x_end=2*pi

N=int(x_end/dx); Nt=int(t_end/dt)
x=zeros(N+2); w0=zeros(N+2); w1=zeros(N+2); w2=zeros(N+2);
w3=zeros(N+2)
w4=zeros(N+2)
#Initial conditions
x[0]=0; w0[0]=0; w1[0]=0; w2[0]=0


for i in range(1,N+1):
    x[i]=x[i-1]+dx
    
    w0[i]=1-cos(x[i])
    w1[i]=w0[i]*(1-r*w0[i+1]+r*w0[i])
    w2[i]=w1[i]*(1-r*w1[i+1]+r*w1[i])
    w3[i]=w2[i]*(1-r*w2[i+1]+r*w2[i])
    w4[i]=w3[i]*(1-r*w3[i+1]+r*w3[i])

figure(figsize=(12,8))
plot(x,w0,x,w1,x,w2,x,w3,x,w4,linewidth=3)
legend((r'$j=0$','$j=1$','$j=2$','$j=3$','$j=4$'),loc='best',fontsize=15)
xlabel(r'$x$',fontsize=20); ylabel(r'$w_{ij}$',fontsize=20)
axis([0,2*pi,0,4])

show()