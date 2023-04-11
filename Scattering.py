# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 14:32:09 2023

@author: Cédric & Jhanene
"""

import numpy as np
import matplotlib.pyplot as plt
from skspatial.objects import Sphere

plt.close('all')


R=200  #Radius of the sphere

pos=np.zeros((1,3))
pos2=np.zeros((1,3))
rightpos=np.zeros((1,3))
dire=np.zeros((1,3))

final_pos=np.zeros((1,3))
final_dir=np.zeros((1,3))



norm = 0
count=0
x2 = np.array([])
y2 = np.array([])
z2 = np.array([])
norm_arr=np.array([])

x2=np.append(x2,0) 
y2=np.append(y2,0) 
z2=np.append(z2,0) 

while (norm<=R) : 

    r=float(np.random.rand(1))
    r1=float(np.random.rand(1))
    r2=float(np.random.rand(1))

    sigma=0.3
    n=0.2

    tau=-np.log(1-r)
    ts=0

    dl=tau/(n*sigma)/100

    phi=2*np.pi*r1
    theta=np.pi*r2

    dire[0,0]=np.sin(theta)*np.cos(phi)
    dire[0,1]=np.sin(theta)*np.sin(phi)
    dire[0,2]=np.cos(theta)
   
    
    while ts<tau :
        
        if norm<=R : 
            for i in range (0,3):
                pos[0,i]=pos[0,i]+dire[0,i]

            norm=np.sqrt(pos[0,0]**2+pos[0,1]**2+pos[0,2]**2)
            ts+=(n*sigma + n*sigma)*dl/2
            norm_arr=np.append(norm_arr,norm)
        else :
            break
        count += 1
    
    if norm<=R : 
        ts2=ts-(n*sigma + n*sigma)*dl/2
        for i in range (0,3):
            pos2[0,i]=pos[0,i]-dire[0,i]
        

        diffts=ts-ts2
        ratio=(tau-ts2)/diffts
    
        for i in range (0,3):
            rightpos[0,i]=pos2[0,i]+ratio*dire[0,i]            
                
        x2=np.append(x2,rightpos[0,0]) 
        y2=np.append(y2,rightpos[0,1]) 
        z2=np.append(z2,rightpos[0,2]) 
    
    else :     
        norm2=norm_arr[count-2]
        for i in range (0,3):
            pos2[0,i]=pos[0,i]-dire[0,i]
        

        diff_norm=norm-norm2
        ratio_norm=(R-norm2)/diff_norm
    
        for i in range (0,3):
            rightpos[0,i]=pos2[0,i]+ratio*dire[0,i]    
        
        x2=np.append(x2,rightpos[0,0])
        y2=np.append(y2,rightpos[0,1]) 
        z2=np.append(z2,rightpos[0,2]) 
        
    final_norm=np.sqrt(rightpos[0,0]**2+rightpos[0,1]**2+rightpos[0,2]**2)

final_pos=rightpos
final_dir=dire
        
fig = plt.figure(figsize=plt.figaspect(1.))

ax = fig.add_subplot(1,1,1 , projection='3d')
ax.grid()

ax.plot(x2,y2,z2)
ax.set_title('1 particle diffusion 3D Plot')
sphere = Sphere([0, 0, 0], R)
sphere.plot_3d(ax, alpha=0.2)
sphere.point.plot_3d(ax, s=100)

# Set axes label
ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=20)
ax.set_zlabel('z', labelpad=20)
plt.show()      
    
            

    
