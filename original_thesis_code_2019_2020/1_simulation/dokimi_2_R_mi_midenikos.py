# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 15:44:59 2019

@author: gd
"""



from __future__ import division 
import math as mth
import numpy as np
import matplotlib.pylab as plt 

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
            azer[:, i] = -z
   return azer  




################################
def R_wfk(w,f,k):
    
    #1
   # w=np.array(3*(np.pi/180))
    #f=np.array(4*(np.pi/180))
    #k=np.array(2*(np.pi/180))

   #px1: #sinw=np.sin(w)
 # px2. : x = np.sin(np.array((0., 30., 45., 60., 90.)) * np.pi / 180. )

       
    #2 edw ayti xrisimopoieitai
    sinw=np.sin(np.array((w))*np.pi/180)
    sinf=np.sin(np.array((f))*np.pi/180)
    sink=np.sin(np.array((k))*np.pi/180)
    
    cosw=np.cos(np.array((w))*np.pi/180)
    cosf=np.cos(np.array((f))*np.pi/180)
    cosk=np.cos(np.array((k))*np.pi/180)
    
    R=np.array([[cosf*cosk,cosw*sink+sinw*sinf*cosk,sinw*sink-cosw*sinf*cosk],
            [-cosf*sink,cosw*cosk-sinw*sinf*sink,sinw*cosk+cosw*sinf*sink],
            [sinf,-sinw*cosf,cosw*cosf]])

    return R
    
    
 
  
##############################
#uvw kai UVW 

def Pimage(c,x0,y0,Pimg):
    Upimg= np.zeros((30,3))
    for i in range(3):
        if i==0:
            Upimg[:, i] =Pimg[:,0]+x0    
        elif i==1:
            Upimg[:, i]=Pimg[:,1]+y0
        else:
            Upimg[:, i] = -c
    return  Upimg  


##########################################################     
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




################################
def ABCDk(x0,y0,z0,P,k,ABCD):
        
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
##################################
    
    
def parallilia(ABCd,k):
    # diereunisi an einai parallila i an taytizontai
     hh=np.zeros((4,4))
     hh[0,0]=ABCd[k,0]/ABCd[(k-1),0]
     hh[0,1]=ABCd[k,1]/ABCd[(k-1),1]
     hh[0,2]=ABCd[k,2]/ABCd[(k-1),2]
     hh[0,3]=ABCd[k,3]/ABCd[(k-1),3]
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

#pinakas strofis:
w=5.0 #moires
f=4.0
k=8.0

R=R_wfk(w,f,k)


#sint/nes toy P (xwroy):
x=15.0 #m
y=25.0
z=4.0

# gia epigeies liyeis : klimaka=-c/(y-y0)    
scale=Y0/(c)
#l=(1/scale)

##########################################################################
            ##kalesmenes sunartiseis:



#1
P=start_points(x,y,z)  

#2
Pimg= Upimage(P,c,X0,Y0,Z0,R) #i sxesi ayti sindeei
# tis snt/nes (xPa,yPa) tis eikonas a me tis epigeies 
#sint/nes XP,YP,ZP twn simeiwn toy P pinaka
#poy proekyye gia ta simeia P toy xwrou. mesw twn stoixeiwn
# eswterikoy pros/smou (x0,y0,c) kai exwterikou
# (X0,Y0,Z0,w,f,k)


#3
Puvw= Pimage(c,x0,y0,Pimg) 
#aytes xrisimopoiw,aytes anaferontai sto sistima eikonosint/nwn!!!!!! 



#xrisimopoiw pinaka Puvw gia na briskomai sto sistima tis
#eikonas ws pros to prwteuon tis.

#4

ABCD=np.zeros((4,4))
k=0
ABCd=ABCDk(X0,Y0,Z0,Puvw,k,ABCD)
#xrisimopoiw pinaka Puvw gia na briskomai sto sistima tis
#eikonas ws pros to prwteuon tis.
k=1
ABCd=ABCDk(X0,Y0,Z0,P,k,ABCD)
k=2
ABCd=ABCDk(0,0,0,P,k,ABCD)
# diereunisi an einai parallila i an taytizontai
k=3
ABCd=ABCDk(0,0,0,Puvw,k,ABCD)
print ABCd
dier=parallilia(ABCd,1)


