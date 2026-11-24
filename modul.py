# -*- coding: utf-8 -*-
"""
Created on Sun May 26 21:47:24 3019

@author: Komiljon
"""
#import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm

from matplotlib.colors import Normalize

#from scipy.interpolate import interp1d

def hello():
    print('Hello, world!')

def fib(n):
    a = b = 1
    for i in range(n - 2):
        a, b = b, a + b
    return b
def grafik(x,c,roa,rop,tm):
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(x[:301],c[:301,  (tm-1)//4],'k',label='$t=450\:c$',lw=1.5)
    ax.plot(x[:301],c[:301,2*(tm-1)//4],'k--',label='$t=900\:c$',lw=2)
    ax.plot(x[:301],c[:301,3*(tm-1)//4],'k-.',label='$t=1350\:c$',lw=2)
    ax.plot(x[:301],c[:301,4*(tm-1)//4],'k:',label='$t=1800\:c$',lw=2)
    #xnew = np.linspace(0, 0.2, 8000)
    #f2 = interp1d(x[:301],c[:301,2700], kind='cubic')

    #ax.plot(xnew,f2(xnew),'k-.',label='$t=1350 c$',lw=1.5)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 1.05, "$c/c_0$",fontsize=14)
    ax.text(ax.get_xlim()[1]*1.01, 0.0, "$x,\: м$", fontsize=14)
    ax.text(ax.get_xlim()[1]/2, 0.98, "$а$",fontsize=14)
    plt.ylim(0.0,1.02)
    #fig.savefig('cc0.eps')
    fig.savefig('cc0.jpg',dpi=600)
    #plt.figure()
    #roa larni chizish
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(x[:301],roa[:301,  (tm-1)//4],'k',label='$t=450\:c$',lw=1.5)
    ax.plot(x[:301],roa[:301,2*(tm-1)//4],'k--',label='$t=900\:c$',lw=2)
    ax.plot(x[:301],roa[:301,3*(tm-1)//4],'k-.',label='$t=1350\:c$',lw=2)
    ax.plot(x[:301],roa[:301,4*(tm-1)//4],'k:',label='$t=1800\:c$',lw=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 0.0512, "$N, \: м^3/м^3$",fontsize=14)
    ax.text(ax.get_xlim()[1]*1.01, 0.0, "$x,\: м$", fontsize=14)
    ax.text(ax.get_xlim()[1]/2, 0.0498, "$б$", fontsize=14)
    plt.ylim(0.0,0.0503)
    #fig.savefig('N.eps')
    fig.savefig('N.jpg',dpi=600)
    #plt.ylim(0.0,0.0103)
    #rop larni chizish
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(x[:301],rop[:301,  (tm-1)//4],'k',label='$t=450\:c$',lw=1.5)
    ax.plot(x[:301],rop[:301,2*(tm-1)//4],'k--',label='$t=900\:c$',lw=2)
    ax.plot(x[:301],rop[:301,3*(tm-1)//4],'k-.',label='$t=1350\:c$',lw=2)
    ax.plot(x[:301],rop[:301,4*(tm-1)//4],'k:',label='$t=1800\:c$',lw=2)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 0.072, "$\u03c1,\: м^3/м^3$",fontsize=14)
    ax.text(ax.get_xlim()[1]*1.01, 0.0, "$x,\:м$", fontsize=14)
    ax.text(ax.get_xlim()[1]/2, 0.0696, "$в$", fontsize=14)
    plt.ylim(0.0,0.0703)
    plt.yticks(0.01*np.arange(0, 8, 1))
    #fig.savefig('ro.eps')
    fig.savefig('ro.jpg',dpi=600)
#c3v1
def grafikc(x,c1,c2,c3):
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)

    ax.plot(x[:301],c2[:301,2700],'k-.',label='$k=0,002$',lw=2)
    ax.plot(x[:301],c1[:301,2700],'k--',label='$k=0,005$',lw=1.5)
    ax.plot(x[:301],c3[:301,2700],'k',label='$k=0,020$',lw=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 1.05, "$c/c_0$",fontsize=14)
    ax.text(ax.get_xlim()[1]*1.01, 0.0, "$x,\: м$", fontsize=14)
    ax.text(ax.get_xlim()[1]/2, 0.98, "$а$",fontsize=14)

    plt.ylim(0.0,1.02)
    #fig.savefig('cc03v1.eps')
    fig.savefig('cc03v1.jpg',dpi=600)

#grafik roa
def grafikroa(x,roa1,roa2,roa3):
    #roa larni chizish
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(x[:301],roa2[:301,2700],'k-.',label='$k=0,002$',lw=2)
    ax.plot(x[:301],roa1[:301,2700],'k--',label='$k=0,005$',lw=1.5)
    ax.plot(x[:301],roa3[:301,2700],'k',label='$k=0,020$',lw=2)    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 0.0512, "$N, \: м^3/м^3$",fontsize=14)
    ax.text(ax.get_xlim()[1]*1.01, 0.0, "$x,\: м$", fontsize=14)
    ax.text(ax.get_xlim()[1]/2, 0.0498, "$б$", fontsize=14)
    plt.ylim(0.0,0.0503)
    #fig.savefig('roa.eps')
    fig.savefig('N3v1.jpg',dpi=600)
    #plt.ylim(0.0,0.0103)
    #rop larni chizish

#rop
def grafikrop(x,rop1,rop2,rop3):
    #rop larni chizish
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(x[:301],rop2[:301,2700],'k-.',label='$k=0,002$',lw=2)
    ax.plot(x[:301],rop1[:301,2700],'k--',label='$k=0,005$',lw=1.5)
    ax.plot(x[:301],rop3[:301,2700],'k',label='$k=0,020$',lw=2)    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 0.075, "$\u03c1_p,  \u0020 м^3/м^3$",fontsize=14)
    ax.text(ax.get_xlim()[1]*1.01, 0.0, "$x,\: м$", fontsize=14)
    ax.text(ax.get_xlim()[1]/2, 0.0696, "$в$", fontsize=14)
    plt.ylim(0.0,0.073)
    plt.yticks(0.01*np.arange(0, 8, 1))
    #fig.savefig('rop.eps')
    fig.savefig('ro3v1.jpg',dpi=600)

##### yangi 3 ta
#c3v1
def grafikcro(x,c1,c2,c3):
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)

    ax.plot(x[:301],c1[:301,2700],'k-.',label='$\u03c1_{1}=0,007$',lw=2)
    ax.plot(x[:301],c2[:301,2700],'k--',label='$\u03c1_{1}=0,010$',lw=1.5)
    ax.plot(x[:301],c3[:301,2700],'k',label='$\u03c1_{1}=0,020$',lw=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 1.05, "$c/c_0$",fontsize=14)
    ax.text(ax.get_xlim()[1]*1.01, 0.0, "$x,\: м$", fontsize=14)
    ax.text(ax.get_xlim()[1]/2, 0.98, "$а$",fontsize=14)

    plt.ylim(0.0,1.02)
    #fig.savefig('cc03v1.eps')
    fig.savefig('cc03v1.jpg',dpi=600)

#grafik roa
def grafikroaro(x,roa1,roa2,roa3):
    #roa larni chizish
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(x[:301],roa1[:301,2700],'k-.',label='$\u03c1_{1}=0,007$',lw=2)
    ax.plot(x[:301],roa2[:301,2700],'k--',label='$\u03c1_{1}=0,010$',lw=1.5)
    ax.plot(x[:301],roa3[:301,2700],'k',label='$\u03c1_{1}=0,020$',lw=2)    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 0.0512, "$N, \: м^3/м^3$",fontsize=14)
    ax.text(ax.get_xlim()[1]*1.01, 0.0, "$x,\: м$", fontsize=14)
    ax.text(ax.get_xlim()[1]/2, 0.0498, "$б$", fontsize=14)
    plt.ylim(0.0,0.0503)
    #fig.savefig('roa.eps')
    fig.savefig('N3v1.jpg',dpi=600)
    #plt.ylim(0.0,0.0103)
    #rop larni chizish

#rop
def grafikropro(x,rop1,rop2,rop3):
    #rop larni chizish
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(x[:301],rop1[:301,2700],'k-.',label='$\u03c1_{1}=0,007$',lw=2)
    ax.plot(x[:301],rop2[:301,2700],'k--',label='$\u03c1_{1}=0,010$',lw=1.5)
    ax.plot(x[:301],rop3[:301,2700],'k',label='$\u03c1_{1}=0,020$',lw=2)    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 0.075, "$\u03c1_p,  \u0020 м^3/м^3$",fontsize=14)
    ax.text(ax.get_xlim()[1]*1.01, 0.0, "$x,\: м$", fontsize=14)
    ax.text(ax.get_xlim()[1]/2, 0.0696, "$в$", fontsize=14)
    plt.ylim(0.0,0.073)
    plt.yticks(0.01*np.arange(0, 8, 1))
    #fig.savefig('rop.eps')
    fig.savefig('ro3v1.jpg',dpi=600)


def grafikroa1(x,roa1,roa2,roa3):
    #roa larni chizish
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(x[:301],roa1[:301,2700],'k',label='$D=0$',lw=1.5)
    ax.plot(x[:301],roa2[:301,2700],'k--',label='$D=2,5\cdot 10^{-6}$',lw=2)
    ax.plot(x[:301],roa3[:301,2700],'k-.',label='$D=5\cdot 10^{-6}$',lw=1.5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 0.0312, "$\u03c1_a, \u0020 м^3/м^3$",fontsize=14)
    ax.text(0.205, 0.0, "$x, м$", fontsize=14)
    ax.text(0.1, 0.0298, "$б$", fontsize=14)
    plt.ylim(0.0,0.0303)
    fig.savefig('roa.eps')
    fig.savefig('roa.jpg')

def grafikroat(x,t,roa1,roa2,roa3):
    #roa larni chizish
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(t,roa3[10,:],'k-.',label='$D=5\cdot 10^{-6}$',lw=2)
    ax.plot(t,roa2[10,:],'k--',label='$D=2,5\cdot 10^{-6}$',lw=2)
    ax.plot(t,roa1[10,:],'k',label='$D=0$',lw=1.5)
    #ax.legend(loc='upper left')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend(loc='lower right')
    ax.text(-0.001, 0.0312, "$\u03c1_a, \u0020 м^3/м^3$",fontsize=14)
    ax.text(1402, 0.0, "$t, c$", fontsize=14)
    #ax.text(0.1, 0.0298, "$б$", fontsize=14)
    plt.ylim(0.0,0.0303)
    fig.savefig('roa.eps')
    fig.savefig('roa.jpg')
    #plt.ylim(0.0,0.0103)
    #rop larni chizish


def grafikrop1(x,rop1,rop2,rop3):
    #rop larni chizish
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(x[:301],rop1[:301,2700],'k',label='$D=0$',lw=1.5)
    ax.plot(x[:301],rop2[:301,2700],'k--',label='$D=2,5\cdot 10^{-6}$',lw=2)
    ax.plot(x[:301],rop3[:301,2700],'k-.',label='$D=5\cdot 10^{-6}$',lw=1.5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 0.075, "$\u03c1_p,  \u0020 м^3/м^3$",fontsize=14)
    ax.text(0.205, 0.0, "$x, м$", fontsize=14)
    ax.text(0.1, 0.0696, "$в$", fontsize=14)
    plt.ylim(0.0,0.073)
    plt.yticks(0.01*np.arange(0, 8, 1))
    fig.savefig('rop.eps')
    fig.savefig('rop.jpg')

def grafikropt(x,t,rop1,rop2,rop3):
    #rop larni chizish
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(t,rop3[10,:],'k-.',label='$D=5\cdot 10^{-6}$',lw=2)
    ax.plot(t,rop2[10,:],'k--',label='$D=2,5\cdot 10^{-6}$',lw=2)
    ax.plot(t,rop1[10,:],'k',label='$D=0$',lw=1.5)
    #ax.legend(loc='upper left')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend(loc='lower right')
    ax.text(-0.001, 0.075, "$\u03c1_p,  \u0020 м^3/м^3$",fontsize=14)
    ax.text(1405, 0.0, "$t, с$", fontsize=14)
    #ax.text(0.1, 0.0696, "$в$", fontsize=14)
    plt.ylim(0.0,0.073)
    plt.yticks(0.01*np.arange(0, 8, 1))
    fig.savefig('ropt.eps')
    fig.savefig('ropt.jpg')



def grafikc1(x,c1,c2,c3):
    fig =plt.figure(figsize=(6, 3))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(x[:301],c1[:301,2700],'k',label='$D=0$',lw=1.5)
    ax.plot(x[:301],c2[:301,2700],'k--',label='$D=2,5\cdot 10^{-6}$',lw=2)
    ax.plot(x[:301],c3[:301,2700],'k-.',label='$D=5\cdot 10^{-6}$',lw=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend()
    ax.text(-0.001, 1.05, "$c/c_0$",fontsize=14)
    ax.text(0.205, 0.0, "$x, м$", fontsize=14)
    ax.text(0.1, 0.98, "$а$",fontsize=14)
    plt.ylim(0.0,1.02)
    fig.savefig('cc0.eps')
    fig.savefig('cc0.jpg')

def grafikct(c,t):
    fig =plt.figure(figsize=(6, 4))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(t,c[8,:],label='0.04',lw=1.5) #,marker="o"
    ax.plot(t,c[40,:],label='$0.2$',lw=1.5)
    ax.plot(t,c[100,:],'k',label='$0.5$',lw=1.5)
    ax.plot(t,c[150,:],label='$0.75$',lw=1.5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    #plt.legend(loc='lower right')
    #ax.text(-0.001, 1.05, "$c/c_0$",fontsize=14)
    #ax.text(1405, 0.0, "$t, c$", fontsize=14)
    #ax.text(0.1, 0.98, "$а$",fontsize=14)
    plt.ylim(0.0,0.012)
    #fig.savefig('cc0.eps')
    fig.savefig('cc0.jpg')

def grafiksigmat(sigma,t):
    fig =plt.figure(figsize=(6, 4))
    plt.style.use('classic')
    ax = plt.subplot(111)
    ax.plot(t,sigma[8,:],label='0.04',lw=1.5)
    ax.plot(t,sigma[40,:],label='$0.2$',lw=1.5)
    ax.plot(t,sigma[100,:],label='$0.5$',lw=1.5)
    ax.plot(t,sigma[150,:],label='$0.75$',lw=1.5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend(loc='center right')
    #ax.text(-0.001, 1.05, "$c/c_0$",fontsize=14)
    #ax.text(1405, 0.0, "$t, c$", fontsize=14)
    #ax.text(0.1, 0.98, "$а$",fontsize=14)
    plt.ylim(0.0,25)
    #fig.savefig('cc0.eps')
    fig.savefig('sigma.jpg')

def graficc3d(c1):
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    w0 = np.arange(0.02,0.41, 0.01) # x - values
    w1 = np.arange(10, 30, 0.5) # y - values
    X, Y = np.meshgrid(w0, w1) # no idea why to do that
    sizeW0 = len(w0)
    sizeW1 = len(w1)
    Z = np.zeros((sizeW1, sizeW0))
    for i in range(0,sizeW0):
        for j in range(0,sizeW1):
            Z[j,i] = c1[i,j]
    minn = Z.min()
    maxx = Z.max()
    #color_map = cm.RdYlGn # reverse the colormap: cm.RdYlGn_r
    #scalarMap = cm.ScalarMappable(norm=Normalize(vmin=min, vmax=max), cmap=color_map)
    #C_colored = scalarMap.to_rgba(Z)
  
    #surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=C_colored, antialiased=True)
  
    #ax.plot_surface(X, Y, Z)
    
    #fig.savefig('FF13D.jpg')
    
    fig1 = plt.figure()
    ax1 = fig1.gca(projection='3d')
    ax1.plot_wireframe(X, Y, Z)
    fig1.savefig('FF3D.jpg')

    #fig = plt.figure(figsize=[15, 10])
    #ax = fig.gca(projection='3d')
    #surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)
    #plt.show()