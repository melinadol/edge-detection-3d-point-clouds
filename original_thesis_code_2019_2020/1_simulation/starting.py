# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\gd\.spyder2\.temp.py
"""


from __future__ import division 
import math as mth
import numpy as np

########
def XYZ0(X0,Y0,Z0):
    XYZ0=np.array([[X0],[Y0],[Z0]])
    return XYZ0
    
    
 #########   
def start_points(x,y,z):
  #επιγειου συστηματος συντνες gia A=pinA(15,25,40)
   azer = np.zeros((30,3))
    
   for i in range(3):
        if i==0:
            azer[:, i] = x
    
        elif i==1: 
            
            for item in range(30): 
                azer[item, i]=y
                y+=1    
       
        elif i==2:
            azer[:, i] = z
   return azer[:,:]    
##################   

 
    
def ABCDk(x0,y0,z0,P,k,AB):
     
    
    
    #Ax+By+Cz+D=0  #eyresi twn A,B,C,D parametrwn mesw pinaka sint/nwn P
    
    a=1
    b=29
    x1=P[a,0]
    x2=P[b,0]
    x3=x0
    
    y1=P[a,1]
    y2=P[b,1]
    y3=y0
    z1=P[a,2]
    z2=P[b,2]
    z3=z0
    A=np.array([[y2-y1,z2-z1],[y3-y1,z3-z1]])
    A=np.linalg.det(A)
   
    B=np.array([[x2-x1,z2-z1],[x3-x1,z3-z1]])
    B=np.linalg.det(B)
     
    C=np.array([[x2-x1,y2-y1],[x3-x1,y3-y1]])
    C=np.linalg.det(C)
    
    D=(-x1*A +y1*B -z1*C )
    

    AB[k,0]=A
    AB[k,1]=B
    AB[k,2]=C
    AB[k,3]=D
    k=k+1
    
    
    return AB
     
    
def ABCD_original_notyet(x,y,z,P):
    pass
    #ABCDp= np.zeros((1000,4))
    #Ax+By+Cz+D=0  #eyresi twn A,B,C,D parametrwn mesw pinaka sint/nwn P
    #i=1
    #a=i
    #b=i+1
    #c=i+2
    #k=1
    #while k<1000:
        
        
    
   
   
   
   
   
#####################     
def Upimage(c,x0,y0,P):
    Upimg= np.zeros((30,3))
    for i in range(3):
        if i==2:
            Upimg[:, i] = -c
        elif i==0:
            Upimg[:, i] =P[:,0]-x0
        else:
            Upimg[:, i]=P[i,1]-y0
    return  Upimg       

######
#uvw kai UVW 
def Pimage(P,c,X0,Y0,R,y0,x0):
    Pimg=np.zeros((30,3))
    for i in range (3):
        if i==0:
            Pimg[:,i]=x0-c*((P[:,0]-X0)*R[0,0]+(P[:,1]-Y0)*R[0,1]+(P[:,2]-Z0)*R[0,2]) / ((P[:,0]-X0)*R[2,0]+(P[:,1]-Y0)*R[2,1]+(P[:,2]-Z0)*R[2,2])
        elif i==1:
            Pimg[:,i]=y0-c*((P[:,0]-X0)*R[1,0]+(P[:,1]-Y0)*R[1,1]+(P[:,2]-Z0)*R[1,2]) /((P[:,0]-X0)*R[2,0]+(P[:,1]-Y0)*R[2,1]+(P[:,2]-Z0)*R[2,2])  
        else :
            Pimg[:,i]=-c
    return Pimg         
 ############################################################ss
  #################################################################
#Stoixeia arxika 
    # global : k, 



P=[]

c=0.24;

X0=100.0
Y0=3.0
Z0=100.0
x0=(36.0/2)*10**-2 ##x0==r
y0=(24.0/2)*10**-2   
z0=-c
    
XYZ0=XYZ0(X0,Y0,Z0)
print XYZ0
#xyz0=XYZ0(x0,y0,z0)

R=np.array([[1,0,0],[0,1,0],[0,0,1]])
x=15.0
y=25.0
z=40.0

######kalesmenes sunartiseis:
#1
P=start_points(x,y,z)  
print P
#2
UPimg= Upimage(c,x0,y0,P)
print UPimg
#3
Pimg= Pimage(P,c,X0,Y0,R,y0,x0) #oi eikonosint/nes toy pinaka P 
#poy proekyye gia ta simeia P toy xwrou.

#4
ABCD= np.zeros((2,4))
print ABCD

k=0
ABCD=ABCDk(x0,y0,z0,P,k,ABCD)
k+=1
ABCD=ABCDk(X0,Y0,Z0,Pimg,k,ABCD)
print ABCD


