# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 15:47:37 2019

@author: gd
"""

from __future__ import division 
import math as mth
import numpy as np
import matplotlib.pylab as plt 


##################################
def ABCDk2(x0,y0,z0,P,k,ABCD):  
        
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
    
###################################
def ABCDk(a,b,X0,Y0,Z0,k,ABCD):
        
    #Ax+By+Cz+D=0  #eyresi twn A,B,C,D parametrwn mesw pinaka sint/nwn P

#simeiwsi        
#2 tyxaia simeia tis eytheias epilegontai, 1,29, kai 3o simeio
#epilegetai i eikona tou simeiou liyis (X0Y0Z0) stin eikona!!!!! 
# 3:(0,0,0) einai oi sintetagmenes tou stin eikona  
 
    
    x1=np.array(a[0,0])
    x2=np.array(b[0,0])
    x3=np.array(X0)
    
    y1=np.array(a[0,1])
    y2=np.array(b[0,1])
    y3=np.array(Y0)
    
    z1=np.array(a[0,2])
    z2=np.array(b[0,2])
    z3=np.array(Z0)
    
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
    

#################### 
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
   
              
                
    
    R=np.dot(R3k,R2f)
    R=np.dot(R,R1w)    
    return R
    
################################
def Upimage(P,c,X0,Y0,Z0,R):   # gia pinakes
    Pimg=np.zeros((30,3))
    for i in range (3):
        if i==0:
            Pimg[:,i]=-c*((P[:,0]-X0)*R[0,0]+(P[:,1]-Y0)*R[0,1]+(P[:,2]-Z0)*R[0,2]) / ((P[:,0]-X0)*R[2,0]+(P[:,1]-Y0)*R[2,1]+(P[:,2]-Z0)*R[2,2])
        elif i==1:
            Pimg[:,i]=-c*((P[:,0]-X0)*R[1,0]+(P[:,1]-Y0)*R[1,1]+(P[:,2]-Z0)*R[1,2]) /((P[:,0]-X0)*R[2,0]+(P[:,1]-Y0)*R[2,1]+(P[:,2]-Z0)*R[2,2])  
        else :
            Pimg[:,i]=-c
    return Pimg 

    
##################################    
def Upimage1(P,c,X0,Y0,Z0,R):     #gia simeia
    Pimg=np.zeros((1,3))
    for i in range (3):
        if i==0:
            Pimg[0,i]=-c*((P[0,0]-X0)*R[0,0]+(P[0,1]-Y0)*R[0,1]+(P[0,2]-Z0)*R[0,2]) / ((P[0,0]-X0)*R[2,0]+(P[0,1]-Y0)*R[2,1]+(P[0,2]-Z0)*R[2,2])
        elif i==1:
            Pimg[0,i]=-c*((P[0,0]-X0)*R[1,0]+(P[0,1]-Y0)*R[1,1]+(P[0,2]-Z0)*R[1,2]) /((P[0,0]-X0)*R[2,0]+(P[0,1]-Y0)*R[2,1]+(P[0,2]-Z0)*R[2,2])  
        else :
            Pimg[0,i]=-c
    return Pimg 
######################    
def start_points(x,y,z):
        # (15,4,25)
#####epigeioy sistimatos sint/nes gia A=pinA(15,25,4)
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
       
###################

    
####################################################################

######kalesmenes sunartiseis:
#1


x=15.0
y=4.0
z=25.0
P=start_points(x,y,z)  
 # gia simeia apodeiji: 
   
A=np.array([[P[0,0],P[0,1],P[0,2]]])
B=np.array([[P[29,0],P[29,1],P[29,2]]])

X0=15.0
Y0=3.0
Z0=10.0
O=np.array([[X0,Y0,Z0]])
c=0.050
R=R_wfk(5.0,4.0,10.0)
pim=Upimage(P,c,X0,Y0,Z0,R)


a=Upimage1(A,c,X0,Y0,Z0,R)
b=Upimage1(B,c,X0,Y0,Z0,R)
ABCd=np.zeros((4,4))

k=0
ABCd=ABCDk(A,B,X0,Y0,Z0,k,ABCd) # A,B ta simeia sto sistima  toy 3D xwroy 

R_1=np.linalg.inv(R)
a_1=np.matrix.transpose(a)
b_1=np.matrix.transpose(b)
xyza_1=np.zeros((3,1))
xyzb_1=np.zeros((3,1))
xyza_1=np.dot(R_1,a_1)
xyzb_1=np.dot(R_1,b_1)
xyza_1=np.matrix.transpose(xyza_1)   # ALLOS TROPOS Y=(xyza_1.T)
xyzb_1=np.matrix.transpose(xyzb_1)
a_0=np.array(O-xyza_1)
b_0=np.array(O-xyzb_1)

k=1
ABCd=ABCDk(a_0,b_0,X0,Y0,Z0,k,ABCd)  #den to xrisimopoiw poythena 
########################################### telos apodeixi simeiwn.

# apodeixi gia pinakes: dimioyrgia yeudosistimatos sint/nwn
pim_1=np.matrix.transpose(pim)
xyzp_1=np.zeros((3,1))
xyzp_1=np.dot(R_1,pim_1)
xyzp_1=(xyzp_1.T)
pim_0=np.array(O-xyzp_1)    #edw egine strofi kai metathesi sint/nwn !!!!

k=2
AK=np.array([[pim_0[0,0],pim_0[0,1],pim_0[0,2]]])
BK=np.array([[pim_0[29,0],pim_0[29,1],pim_0[29,2]]])
ABCd=ABCDk(AK,BK,X0,Y0,Z0,k,ABCd)


#  PINAKAS ME YEUDOSISTIMA TELIKOS === pim_0



hh=np.zeros((4,1))

hh[0,0]=np.array(ABCd[0,0]/ABCd[2,0])
#hh[1,0]=np.array(ABCd[0,1]/ABCd[2,1])
#hh[2,0]=np.array(ABCd[0,2]/ABCd[2,2])
hh[3,0]=np.array(ABCd[0,3]/ABCd[2,3])

#################################################################################################
###############################################################
##############################################################


def elegxos_simeiwn_se_epipeda(ABCd,UVW,k):
    debug=np.zeros((3,30))
    for i in range(30):
        debug[k,i]=ABCd[k,0]*UVW[i,0]+ABCd[k,1]*UVW[i,1]+ABCd[k,2]*UVW[i,2]+ABCd[k,3]
    return debug
    
    
for k in range(3):
    if k==0:
        debug=elegxos_simeiwn_se_epipeda(ABCd,P,k)
    elif k==1:
        debug=elegxos_simeiwn_se_epipeda(ABCd,pim_0,k)
    elif k==2:
        debug=elegxos_simeiwn_se_epipeda(ABCd,pim_0,k)
        
     
###  APODEIXTHIKE OTI OLA TA SIMEIA EPALITHEVOUN TA EPIPEDA XWRIS SFALMA 




pin_prosd=np.zeros((29,4))
k=0
for i in range(29):
    a=np.array([[pim_0[i,0],pim_0[i,1],pim_0[i,2]]])
    b=np.array([[pim_0[i+1,0],pim_0[i+1,1],pim_0[i+1,2]]])
    pin_prosd=ABCDk(a,b,X0,Y0,Z0,k,pin_prosd)
    k+=1

       
        