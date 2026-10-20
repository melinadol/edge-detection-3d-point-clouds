# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 16:08:45 2019

@author: gd
"""
from __future__ import division 
#import math as mth
import numpy as np
#import matplotlib.pylab as plt 
def sfalma_sint_nwn_kai_epipeda(P):
    pass
    


def elegxos_simeiwn_se_epipeda(ABCd,UVW,k):
    debug=np.zeros(30)
    for i in range(30):
        debug[i]=ABCd[k,0]*UVW[i,0]+ABCd[k,1]*UVW[i,1]+ABCd[k,2]*UVW[i,2]+ABCd[k,3]
    return debug
    
    
def ABCD_every_point(X0,Y0,Z0,Pimg,m,ABCDp):
    # DEN  EXEI XRISIMOPOIITHEI
    x3=np.array(Pimg[m,0])
    x2=np.array(0)
    x1=np.array(X0)
    
    y3=np.array(Pimg[m,1])
    y2=np.array(0)
    y1=np.array(Y0)
    
    z3=np.array(Pimg[m,2])
    z2=np.array(0)
    z1=np.array(Z0)
    
    A=np.array([[(y2-y1),(z2-z1)],[(y3-y1),(z3-z1)]])
   
    A=np.linalg.det(A)
    
    
   
    B=np.array([[(x2-x1),(z2-z1)],[(x3-x1),(z3-z1)]])
    
    B=np.linalg.det(B)
    
    
       
    C=np.array([[(x2-x1),(y2-y1)],[(x3-x1),(y3-y1)]])
    
    C=np.linalg.det(C)
    
    
    
    D=(-x1*A +y1*B -z1*C )
    

    ABCDp[m,0]=A
    ABCDp[m,1]=-B
    ABCDp[m,2]=C
    ABCDp[m,3]=D
    
    return ABCDp
        
    
    
#####################################################
def ABCDk(X0,Y0,Z0,P,k,ABCD):
        
    #Ax+By+Cz+D=0  #eyresi twn A,B,C,D parametrwn mesw pinaka sint/nwn P

#simeiwsi        
#2 tyxaia simeia tis eytheias epilegontai, 1,29, kai 3o simeio
#epilegetai i eikona tou simeiou liyis (X0Y0Z0) stin eikona!!!!! 
# 3:(0,0,0) einai oi sintetagmenes tou stin eikona  
 
    a=0
    b=29
    x1=np.array(P[a,0])
    x2=np.array(P[b,0])
    x3=np.array(x0)
    
    y1=np.array(P[a,1])
    y2=np.array(P[b,1])
    y3=np.array(y0)
    
    z1=np.array(P[a,2])
    z2=np.array(P[b,2])
    z3=np.array(z0)
    
    A=np.array([[(y2-y1),(z2-z1)],[(y3-y1),(z3-z1)]])
    
    A=np.linalg.det(A)
    
    
   
    B=np.array([[(x2-x1),(z2-z1)],[(x3-x1),(z3-z1)]])
    
    B=np.linalg.det(B)
     
    
    
    C=np.array([[(x2-x1),(y2-y1)],[(x3-x1),(y3-y1)]])
    
    C=np.linalg.det(C)
    
     
    D=(-x1*A +y1*B -z1*C )
    

    ABCD[k,0]=A
    ABCD[k,1]=-B
    ABCD[k,2]=C
    ABCD[k,3]=D
    
    deb=A*x2-B*y2+C*z2+D
    print "debug is"
    print deb
    
    return ABCD[:,:]
    

    
    
    
    
    
    
############################################################################
def matmul(A,B):
    A=np.array(A)
    B=np.array(B)
    if (len(A[0]))==(len(B)):
        C=np.zeros((len(A),len(B[0])))
        # iterate through rows of A
        for i in range(len(A)):
               # iterate through columns of B
               for j in range(len(B[0])):
                   # iterate through rows of B
                   for k in range(len(B)):
                       C[i,j] += A[i,k] * B[k,j]
        
        
    else:
         C='den ginetai pollaplasiasmos xxxxxxx'
         
        
    return C
    
##############

def start_points(x,y,z):
        # (15,4,25)
  #επιγειου συστηματος συντ/νες gia A=pinA(15,25,4)
   azer = np.zeros((30,3))
    
   for i in range(3):
        if i==0:
            azer[:, i] = x
            
        elif i==1:
            azer[:, i] = y
            
        elif i==2: 
            for item in range(30): 
                #ousiastika ta z anaferontai sta diaforetika yyi sto xwro
                azer[item, i]=z
                z=z+0.3    
            
       
        
   return azer
   
###############################################   
   
def R_wfk(w,f,k):
    
    #1
   # w=np.array(3*(np.pi/180))
    #f=np.array(4*(np.pi/180))
    #k=np.array(2*(np.pi/180))

   #px1: #sinw=np.sin(w)
 # px2. : x = np.sin(np.array((0., 30., 45., 60., 90.)) * np.pi / 180. )

       
   #2 edw ayti xrisimopoieitai
    
    
    sinw=np.sin(np.array((w))*np.pi/180)
    cosw=np.cos(np.array((w))*np.pi/180)
    
    sinf=np.sin(np.array((f))*np.pi/180)
    cosf=np.cos(np.array((f))*np.pi/180)
    
    
    sink=np.sin(np.array((k))*np.pi/180)
    cosk=np.cos(np.array((k))*np.pi/180)   
   
    R1w=np.array([[1,0,0],
                  [0,cosw,sinw],
                  [0,-sinw,cosw]])
                  
 
   
    R2f=np.array([[cosf,0,-sinf],
                  [0,1,0],
                  [sinf,0,cosf]])
   
   
    R3k=np.array([[cosk,sink,0],
                  [-sink,cosk,0],
                  [0,0,1]])
   
              
                
    
    R=matmul(R3k,R2f)
    R=matmul(R,R1w)    
    return R
        
    
def strofi_metath_UVW(P,X0,Y0,Z0,Rst): #KAI PARALLILI METATHESI MAZI (X-X0,Y-Y0,Z-Z0)
    P_XyZ0=np.zeros((30,3))
    P_XyZ0=np.array(P-[X0,Y0,Z0])
    print P_XyZ0
    UVW=np.zeros((30,3))
    for i in range(30):
        UVW[i,0]=Rst[0,0]*P_XyZ0[i,0]+Rst[0,1]*P_XyZ0[i,1]+Rst[0,2]*P_XyZ0[i,2]
        UVW[i,1]=Rst[1,0]*P_XyZ0[i,0]+Rst[1,1]*P_XyZ0[i,1]+Rst[1,2]*P_XyZ0[i,2]
        UVW[i,2]=Rst[2,0]*P_XyZ0[i,0]+Rst[2,1]*P_XyZ0[i,1]+Rst[2,2]*P_XyZ0[i,2]
    return UVW    

    
def strofi_xwris_metath_UVW(P_XyZ0,Rst):  # DEN XRISISMOPOIEITAI
    UVW=np.zeros((30,3))
    for i in range(30):
        UVW[i,0]=Rst[0,0]*P_XyZ0[i,0]+Rst[0,1]*P_XyZ0[i,1]+Rst[0,2]*P_XyZ0[i,2]
        UVW[i,1]=Rst[1,0]*P_XyZ0[i,0]+Rst[1,1]*P_XyZ0[i,1]+Rst[1,2]*P_XyZ0[i,2]
        UVW[i,2]=Rst[2,0]*P_XyZ0[i,0]+Rst[2,1]*P_XyZ0[i,1]+Rst[2,2]*P_XyZ0[i,2]
    return UVW        
    
    
def Upimage(P,c,X0,Y0,Z0,R):
    Pimg=np.zeros((30,3))
    for i in range (3):
        if i==0:
            Pimg[:,i]=-c*((P[:,0]-X0)*R[0,0]+(P[:,1]-Y0)*R[0,1]+(P[:,2]-Z0)*R[0,2]) / ((P[:,0]-X0)*R[2,0]+(P[:,1]-Y0)*R[2,1]+(P[:,2]-Z0)*R[2,2])
        elif i==1:
            Pimg[:,i]=-c*((P[:,0]-X0)*R[1,0]+(P[:,1]-Y0)*R[1,1]+(P[:,2]-Z0)*R[1,2]) /((P[:,0]-X0)*R[2,0]+(P[:,1]-Y0)*R[2,1]+(P[:,2]-Z0)*R[2,2])  
        else :
            Pimg[:,i]=-c
    return Pimg 



   
def klimaka(S1,S2):
    pass

#
def Pimage(c,x0,y0,Pimg):
    Upimg= np.zeros((30,3))
    for i in range(3):
        if i==2:
            Upimg[:, i] = -c
        elif i==0:
            Upimg[:, i] =Pimg[:,0]+x0
        else:
            Upimg[:, i]=Pimg[:,1]+y0
    return  Upimg          
  
   
x=15.0
y=4.0
z=25.0


######kalesmenes sunartiseis:
#1
P=start_points(x,y,z)  


w=90.0 #moires STROFI  SYSTIMATWN SINTETAGMENWN
f=0.0
k=0.0
#strofes
Rstrofi=R_wfk(4.,5.,8.0)

X0=20.0
Y0=3.0
Z0=10.0

UVW_p=strofi_metath_UVW(P,X0,Y0,Z0,Rstrofi)
XYZ_p=strofi_xwris_metath_UVW(P,Rstrofi)

#     EFARMOGI SIGGRAMIKOTHTAS
c=0.24*10**-1
Rwfk=R_wfk(5.0,4.0,8.0)   #moires  ANAFERETAI STI STROFI GIA SIGGRAMIKOTHTA
#      w=5,  f=4  ,  k=8.0

x0=0
y0=0
z0=0
uvw_P=Upimage(P,c,X0,Y0,Z0,Rwfk)
uvw_P2=Upimage(UVW_p,c,X0,Y0,Z0,Rwfk)
Pimg=Pimage(c,x0,y0,uvw_P)

ABCd=np.zeros((6,4))

k=0
ABCd=ABCDk(X0,Y0,Z0,UVW_p,k,ABCd)
k=1
ABCd=ABCDk(0,0,0,UVW_p,k,ABCd)
k=2
ABCd=ABCDk(X0,Y0,Z0,uvw_P,k,ABCd)
k=3
ABCd=ABCDk(0,0,0,uvw_P,k,ABCd)
k=4
ABCd=ABCDk(X0,Y0,Z0,uvw_P2,k,ABCd)

ABCDevU=np.zeros((30,4))
ABCDevu=np.zeros((30,4))
for m in range(30):
    ABCDevU=ABCD_every_point(X0,Y0,Z0,UVW_p,m,ABCDevU)
    ABCDevu=ABCD_every_point(X0,Y0,Z0,uvw_P,m,ABCDevu)

############################################################################
for k in range(4):
    if k==0:
        debug0=elegxos_simeiwn_se_epipeda(ABCd,UVW_p,k)
    elif k==1:
        debug1=elegxos_simeiwn_se_epipeda(ABCd,UVW_p,k)
    elif k==2:
        debug2=elegxos_simeiwn_se_epipeda(ABCd,uvw_P,k)
    else:
        debug3=elegxos_simeiwn_se_epipeda(ABCd,uvw_P,k)
        
###  APODEIXTHIKE OTI OLA TA SIMEIA EPALITHEVOUN TA EPIPEDA XWRIS SFALMA 
randomS=np.linspace(-0.01,0.01,90)  
UVW_fake1=np.zeros((30,3)) 
UVW_fake2=np.zeros((30,3)) 
for j in range(30): 
    UVW_fake1[j,0]=UVW_p[j,0]+randomS[j]
    UVW_fake2[j,0]=UVW_p[j,0]-randomS[j]
for j in range(30):
    UVW_fake1[j,1]=UVW_p[j,1]+randomS[j+30]
    UVW_fake2[j,1]=UVW_p[j,1]-randomS[j+30]   
for j in range(30):
        UVW_fake1[j,2]=UVW_p[j,2]+randomS[j+60]
        UVW_fake2[j,2]=UVW_p[j,2]-randomS[j+60]
       
    
        
        
    


epiped=np.zeros((4,4))
k=0
epiped=ABCDk(X0,Y0,Z0,UVW_fake1,k,epiped)
#epiped=epiped-ABCd[k,:]
k=1
epiped=ABCDk(0,0,0,UVW_fake1,k,epiped)  
k=2
epiped=ABCDk(X0,Y0,Z0,UVW_fake2,k,epiped)
#epiped=epiped-ABCd[k,:]
k=3
epiped=ABCDk(0,0,0,UVW_fake2,k,epiped)  

diaf_epip=np.zeros((4,4))
for k in range(4):
    if k<2:
        diaf_epip[k,:]=ABCd[k,:]-epiped[k,:]  
    else:
        diaf_epip[k,:]=ABCd[k-2,:]-epiped[k,:]
    
        