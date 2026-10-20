# -*- coding: utf-8 -*-
"""
Created on Thu Apr 04 01:30:24 2019

@author: gd
"""
from __future__ import division 
#import math as mth
import numpy as np
#import matplotlib.pylab as plt
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
    
    
    
def Upimage(P,c,X0,Y0,Z0,R):
    Pimg=np.zeros((1,3))
    for i in range (3):
        if i==0:
            Pimg[0,i]=-c*((P[0,0]-X0)*R[0,0]+(P[0,1]-Y0)*R[0,1]+(P[0,2]-Z0)*R[0,2]) / ((P[0,0]-X0)*R[2,0]+(P[0,1]-Y0)*R[2,1]+(P[0,2]-Z0)*R[2,2])
        elif i==1:
            Pimg[0,i]=-c*((P[0,0]-X0)*R[1,0]+(P[0,1]-Y0)*R[1,1]+(P[0,2]-Z0)*R[1,2]) /((P[0,0]-X0)*R[2,0]+(P[0,1]-Y0)*R[2,1]+(P[0,2]-Z0)*R[2,2])  
        else :
            Pimg[0,i]=-c
    return Pimg 


    
    
A=np.array([[10.0,10.0,10.0]])
B=np.array([[15.0,25.0,5.0]])
X0=15.0
Y0=15.0
Z0=30.0
O=np.array([[X0,Y0,Z0]])
c=0.050
R=R_wfk(5.0,4.0,10.0)
a=Upimage(A,c,X0,Y0,Z0,R)
b=Upimage(B,c,X0,Y0,Z0,R)
ABCd=np.zeros((4,4))
k=0
ABCd=ABCDk(A,B,X0,Y0,Z0,k,ABCd)
k=1
#ABCd=ABCDk(a,b,X0,Y0,Z0,k,ABCd)

R_1=np.linalg.inv(R)
a_1=np.matrix.transpose(a)
b_1=np.matrix.transpose(b)
xyz_1=np.zeros((3,1))
xyza_1=np.dot(R_1,a_1)
xyzb_1=np.dot(R_1,b_1)
xyza_1=np.matrix.transpose(xyza_1)   # ALLOS TROPOS Y=(xyza_1.T)
xyzb_1=np.matrix.transpose(xyzb_1)


a_0=np.array(O-xyza_1)
b_0=np.array(O-xyzb_1)

k=2
ABCd=ABCDk(a_0,b_0,X0,Y0,Z0,k,ABCd)


hh=np.zeros(4)
hh[0]=ABCd[0,0]/ABCd[2,0]
hh[1]=ABCd[0,1]/ABCd[2,1]
hh[2]=ABCd[0,2]/ABCd[2,2]
hh[3]=ABCd[0,3]/ABCd[2,3]



