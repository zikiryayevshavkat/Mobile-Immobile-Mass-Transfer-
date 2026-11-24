# -*- coding: utf-8 -*-
"""
Created on Thu May 23 18:26:55 2019

@author: Komiljon
"""
import numpy as np
#import matplotlib.pyplot as plt
import copy
import modul
#boshlang'ich ma'lumotlar
#beta=0.005

m0=0.3
m2=1#0.15
v=0.0001
c00=0.05

D=2.5e-6

r0=0.1
ra0=0.05
rp0=0.07

rp1=0#.008
ra1=0.001

k1=0.005#*0.4*100
k2=0.005#*0.4*100

betar=0.003
betaa=0.003*0.7/0.3*1.11
betad=0.005

nn=301
tm=3601 

#2700
#N=range(nn)
#print(N[0])
tmax=range(1,tm) #vaqt
tau=0.5
h=0.001
#Procedure Nach_Dan;{Boshlang'ich va chegaraviy shartlar}
N=np.zeros((nn,tm), dtype=float)
ro=np.zeros((nn,tm), dtype=float)
c=np.zeros((nn,tm), dtype=float)
alpha=np.zeros(nn, dtype=float)
bet=np.zeros(nn, dtype=float)

c[0,:]=c00
#  /// bettaRo ni fisoblash funksiyasi

def bettaRo(ro,cij):
   if  ro<=rp1: 
     return betar*cij
   elif ro<rp0:
     return betaa*cij-betad*ro
   else: 
     return 0.0
#  /// bettaN ni fisoblash funksiyasi
def bettaN(N,cij):
   #if  ra<=ra1: 
    # return betar*v*cij
   #elif ra<ra0:
     return k1*cij-k2*N
   #else: 
     #return 0.0
def progonka():
    A=tau*D/(h*h)+tau*v/h
    B=2*tau*D/(h*h)+tau*v/h+m0
    E=tau*D/(h*h)
    alpha[1]=0
    bet[1]=c00
    for i in range(1,nn-1):
        F=m0*c[i,j-1]-(ro[i,j]-ro[i,j-1])-m2*(N[i,j]-N[i,j-1])
        bet[i+1]=(F+A*bet[i])/(B-A*alpha[i])
        alpha[i+1]=E/(B-A*alpha[i])
    c[nn-1,j]=(F+A*bet[nn-1])/(B-A*alpha[nn-1])
    for i in reversed(range(nn-1)):
        c[i,j]=alpha[i+1]*c[i+1,j]+bet[i+1]

#Boshlandi
for j in tmax:
     for i in range(nn):
         N[i,j]=N[i,j-1]+tau*bettaN(N[i,j-1],c[i,j-1])
         ro[i,j]=ro[i,j-1]+tau*bettaRo(ro[i,j-1],c[i,j-1])
     progonka()
        #for i in range(1,nn):
         #c[i,j]=(m0*h*c[i,j-1]+v*tau*c[i-1,j]-h*(roa[i,j]-roa[i,j-1]+rop[i,j]-rop[i,j-1]))/(m0*h+v*tau)
#grafiklar
c=c/c00

#modul.grafik(x,c,roa,rop,tm)
roa1=copy.copy(N)
rop1=copy.copy(ro)
c1=copy.copy(c)

## 2- qiymat

m0=0.3
m2=1#0.15
v=0.0001
c00=0.05

D=2.5e-6

r0=0.1
ra0=0.05
rp0=0.07

rp1=0#.008
ra1=0.001

k1=0.005*0.4#*100
k2=0.005*0.4#*100

betar=0.003
betaa=0.003*0.7/0.3*1.11
betad=0.005


#Boshlandi

roa=np.zeros((nn,tm), dtype=float)
rop=np.zeros((nn,tm), dtype=float)
c=np.zeros((nn,tm), dtype=float)
c[0,:]=c00

for j in tmax:
     for i in range(nn):
         N[i,j]=N[i,j-1]+tau*bettaN(N[i,j-1],c[i,j-1])
         ro[i,j]=ro[i,j-1]+tau*bettaRo(ro[i,j-1],c[i,j-1])
     progonka()
#grafiklar
c=c/c00
#modul.grafik(x,c,roa,rop,tm)
roa2=copy.copy(N)
rop2=copy.copy(ro)
c2=copy.copy(c)
## 3- qiymat

m0=0.3
m2=1#0.15
v=0.0001
c00=0.05

D=2.5e-6

r0=0.1
ra0=0.05
rp0=0.07

rp1=0#.008
ra1=0.001

k1=0.005*0.4*100
k2=0.005*0.4*100

betar=0.003
betaa=0.003*0.7/0.3*1.11
betad=0.005



#Boshlandi

roa=np.zeros((nn,tm), dtype=float)
rop=np.zeros((nn,tm), dtype=float)
c=np.zeros((nn,tm), dtype=float)
c[0,:]=c00

for j in tmax:
     for i in range(nn):
         N[i,j]=N[i,j-1]+tau*bettaN(N[i,j-1],c[i,j-1])
         ro[i,j]=ro[i,j-1]+tau*bettaRo(ro[i,j-1],c[i,j-1])
     progonka()
#grafiklar
c=c/c00
#modul.grafik(x,c,roa,rop,tm)
roa3=copy.copy(N)
rop3=copy.copy(ro)
c3=copy.copy(c)
## 2- qiymat
x = np.linspace(0, h*(nn-1), num=nn, endpoint=True)
modul.grafikc(x,c1,c2,c3)
modul.grafikroa(x,roa1,roa2,roa3)
modul.grafikrop(x,rop1,rop2,rop3)
#tt = np.linspace(0, tau*(tm-1), num=tm, endpoint=True)
#modul.grafikropt(0.01,tt,rop1,rop2,rop3)
#modul.grafikroat(0.01,tt,roa1,roa2,roa3)

