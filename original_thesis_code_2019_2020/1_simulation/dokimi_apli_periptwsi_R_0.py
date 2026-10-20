# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 15:34:39 2018

@author: gd
"""

from __future__ import division 
#import math as mth
import numpy as np

########
def XYZ0(X0,Y0,Z0):
    XYZ0=np.array([[X0],[Y0],[Z0]])
    return XYZ0[:,:]
    
    
 #########   
def start_points(x,y,z):
  #επιγειου συστηματος συντ/νες gia A=pinA(15,25,40)
   azer = np.zeros((30,3))
    
   for i in range(3):
        if i==0:
            azer[:, i] = x
    
        elif i==1: 
            
            for item in range(30): 
                #ousiastika ta y anaferontai sta diaforetika yyi sto xwro
                azer[item, i]=y
                y=y+0.3    
       
        elif i==2:
            azer[:, i] = z
   return azer[:,:]    
##################   


 

    
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
            Upimg[:, i] =Pimg[:,0]-x0
        else:
            Upimg[:, i]=Pimg[:,1]-y0
    return  Upimg          

######
#uvw kai UVW 
def Pimage(P,c,X0,Y0,Z0,R,y0,x0):
    Pimg=np.zeros((30,3))
    for i in range (3):
        if i==0:
            Pimg[:,i]=x0-c*((P[:,0]-X0)*R[0,0]+(P[:,1]-Y0)*R[0,1]+(P[:,2]-Z0)*R[0,2]) / ((P[:,0]-X0)*R[2,0]+(P[:,1]-Y0)*R[2,1]+(P[:,2]-Z0)*R[2,2])
        elif i==1:
            Pimg[:,i]=y0-c*((P[:,0]-X0)*R[1,0]+(P[:,1]-Y0)*R[1,1]+(P[:,2]-Z0)*R[1,2]) /((P[:,0]-X0)*R[2,0]+(P[:,1]-Y0)*R[2,1]+(P[:,2]-Z0)*R[2,2])  
        else :
            Pimg[:,i]=-c
    return Pimg[:,:]     




def ABCDk(x0,y0,z0,L,k,ABCD):
        
    #Ax+By+Cz+D=0  #eyresi twn A,B,C,D parametrwn mesw pinaka sint/nwn P

#simeiwsi        
#2 tyxaia simeia tis eytheias epilegontai, 1,29, kai 3o simeio
#epilegetai i eikona tou simeiou liyis (X0Y0Z0) stin eikona!!!!! 
# 3:(0,0,0) einai oi sintetagmenes tou stin eikona  
    
    a=0
    b=29
    x1=np.array(L[a,0])
    x2=np.array(L[b,0])
    x3=np.array(x0)
    
    y1=np.array(L[a,1])
    y2=np.array(L[b,1])
    y3=np.array(y0)
    
    z1=np.array(L[a,2])
    z2=np.array(L[b,2])
    z3=np.array(z0)
    
    A=np.array([[y2-y1,z2-z1],[y3-y1,z3-z1]])
    A=np.linalg.det(A)
   
    B=np.array([[x2-x1,z2-z1],[x3-x1,z3-z1]])
    B=np.linalg.det(B)
     
    C=np.array([[x2-x1,y2-y1],[x3-x1,y3-y1]])
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
########################################    
def parallilia(ABCD):
    # diereunisi an einai parallila i an taytizontai
     hh=np.zeros(4)
     hh[0]=ABCd[1,0]/ABCd[0,0]
     hh[1]=ABCd[1,1]/ABCd[0,1]
     hh[2]=ABCd[1,2]/ABCd[0,2]
     hh[3]=ABCd[1,3]/ABCd[0,3]
     return hh
    
    
    
    
         
 ############################################################ss
  #################################################################


#Stoixeia arxika 
    # global : k, 



P=[]

c=0.24*10**-1

X0=20.0
Y0=3.0
Z0=10.0

x0=0
y0=0
z0=-c

XYZ0=XYZ0(X0,Y0,Z0)
print XYZ0
#xyz0=XYZ0(x0,y0,z0)

R=np.eye(3)
x=15.0
y=25.0
z=4.0

######kalesmenes sunartiseis:
#1
P=start_points(x,y,z)  
print P
#2

Pimg=Pimage(P,c,X0,Y0,Z0,R,y0,x0) 
#i sxesi ayti sindeei
# tis snt/nes (xPa,yPa) tis eikonas a me tis epigeies 
#sint/nes XP,YP,ZP twn simeiwn toy P pinaka
#poy proekyye gia ta simeia P toy xwrou. mesw twn stoixeiwn
# eswterikoy pros/smou (x0,y0,c) kai exwterikou
# (X0,Y0,Z0,w,f,k)


#3
Puvw= Upimage(c,x0,y0,Pimg)
print Puvw

#xrisimopoiw pinaka Puvw gia na briskomai sto sistima tis
#eikonas ws pros to prwteuon tis.
#4

ABCd=np.zeros((4,4))
k=0
ABCd=ABCDk(X0,Y0,Z0,Puvw,k,ABCd)
k=1
ABCd=ABCDk(X0,Y0,Z0,P,k,ABCd)
k=2
ABCd=ABCDk(0,0,0,P,k,ABCd)
k=3
ABCd=ABCDk(0,0,0,Puvw,k,ABCd)

print ABCd



