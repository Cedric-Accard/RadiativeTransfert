# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 14:10:07 2023

@author: Cédric & Jhanene
"""


import numpy as np
import matplotlib.pyplot as plt
from skspatial.objects import Sphere


photons=5 #Numbers of photons
R=200  #Radius of the sphere

origin=[1e-200,1e-200,1e-200] #Source Origin 

#Initialization of arrays for telescope projection 
z_image=np.zeros(photons)
x_image=np.zeros(photons)

#Initialization of arrays for keeping track of all partciles positions
x_trajectories=np.zeros(shape=(photons,100))
y_trajectories=np.zeros(shape=(photons,100))
z_trajectories=np.zeros(shape=(photons,100))



for j in range(photons) : #For each particle do the scattering loop
    
    #Initialization of arrays for position and direction
    pos=np.zeros((1,3))
    for k in range(0,3):
        pos[0,k]=origin[k]
    pos2=np.zeros((1,3))
    rightpos=np.zeros((1,3))
    dire=np.zeros((1,3))

    final_pos=np.zeros((1,3))
    final_dir=np.zeros((1,3))
    
    norm = 0
    count = 0
    #Initialization of arrays to store future values of positions and add start with the source origin
    x2 = np.array([])
    y2 = np.array([])
    z2 = np.array([])
    norm_arr=np.array([])
    x2=np.append(x2,origin[0]) 
    y2=np.append(y2,origin[1]) 
    z2=np.append(z2,origin[2]) 
    
    while (norm<=R) :  #While we are inside the oject : do the scattering
    
        #Random number generators for the scattering of each particles
        r=float(np.random.rand(1))
        r1=float(np.random.rand(1))
        r2=float(np.random.rand(1))
        
        sigma=0.3 #Cross-Section
        n=0.2 #Density
    
        tau=-np.log(1-r) #Optical Depth
        ts=0
        dl=tau/(n*sigma)/100 #Integration path
        
        # Definition of the direction vector    
        phi=2*np.pi*r1
        theta=np.pi*r2
        dire[0,0]=np.sin(theta)*np.cos(phi)
        dire[0,1]=np.sin(theta)*np.sin(phi)
        dire[0,2]=np.cos(theta)
       
        while ts<tau :
            
            if norm<=R :                        #Going forward by adding dire vector to the postion while inside the object
                for i in range (0,3):
                    pos[0,i]=pos[0,i]+dire[0,i]
    
                norm=np.sqrt((pos[0,0]-origin[0])**2+(pos[0,1]-origin[1])**2+(pos[0,2]-origin[2])**2)
                ts=ts+(n*sigma + n*sigma)*dl/2
                norm_arr=np.append(norm_arr,norm)
            else :
                break
            count += 1
        
        if norm<=R : 
            ts2=ts-(n*sigma + n*sigma)*dl/2     #Correction of the overshoot
            for i in range (0,3):               #
                pos2[0,i]=pos[0,i]-dire[0,i]    #
            
    
            diffts=ts-ts2
            ratio=(tau-ts2)/diffts
        
            for i in range (0,3):
                rightpos[0,i]=pos2[0,i]+ratio*dire[0,i] #Correction to the position to get forward as close as possible to the limit of the object             
                    
            x2=np.append(x2,rightpos[0,0]) # Adding the position to a list to keep track of them
            y2=np.append(y2,rightpos[0,1]) #
            z2=np.append(z2,rightpos[0,2]) #
        
        else :     
            norm2=norm_arr[count-2]
            for i in range (0,3):
                pos2[0,i]=pos[0,i]-dire[0,i]    
            
    
            diff_norm=norm-norm2            #Calculation to get the ratio of dire vector to add to the position to correct it
            ratio_norm=(R-norm2)/diff_norm
        
            for i in range (0,3):
                rightpos[0,i]=pos2[0,i]+ratio_norm*dire[0,i]    #Correction to the position to get backward as close as possible to the limit of the object
            
            x2=np.append(x2,rightpos[0,0]) # Adding the position to a list to keep track of them
            y2=np.append(y2,rightpos[0,1]) # 
            z2=np.append(z2,rightpos[0,2]) #
            
        final_norm=np.sqrt(rightpos[0,0]**2+rightpos[0,1]**2+rightpos[0,2]**2) #Last norm 
    
    final_pos=rightpos # Last position of the particle at the limit of the object
    final_dir=dire # Last direction of the particle at the limit of the object    
    
    theta=np.arccos(final_pos[0,2]/final_norm) ##Definition of the theta angle
    phi=np.arctan2(final_pos[0,1], final_pos[0,0])        

    x_image[j]=final_norm*np.cos(phi)*np.sin(theta) #Definition of the horizontal offset
    z_image[j]=final_pos[0,2] #final_norm*np.sin(phi)*np.sin(theta) #Z projection of the particle postion     

    
    for k in range(len(x2)):
        x_trajectories[j,k]=x2[k]
    for k in range(len(y2)):
        y_trajectories[j,k]=y2[k]
    for k in range(len(z2)):
        z_trajectories[j,k]=z2[k]

fig = plt.figure(figsize=plt.figaspect(2.))

#Ploting the sphere and the trajectories
ax = fig.add_subplot(1,2,1 , projection='3d')
ax.grid()

#Remove useless points
x_trajectories[ x_trajectories==0 ] = np.nan
y_trajectories[ y_trajectories==0 ] = np.nan
z_trajectories[ z_trajectories==0 ] = np.nan

#Trajectories
for l in range(0,photons) : 
    ax.plot(x_trajectories[l,],y_trajectories[l,],z_trajectories[l,])
    
ax.set_title('many particle diffusion 3D Plot')
#Sphere
sphere = Sphere([origin[0], origin[1], origin[2]], 200)
sphere.plot_3d(ax, alpha=0.2)
sphere.point.plot_3d(ax, s=100)

# Set axes label
ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=20)
ax.set_zlabel('z', labelpad=20)

#Telescope
ax = fig.add_subplot(1, 2, 2)
ax.scatter(x_image,z_image,  marker=".")

yabs_max = abs(max(ax.get_ylim(), key=abs))
ax.set_ylim(ymin=-yabs_max, ymax=yabs_max)

xabs_max = abs(max(ax.get_xlim(), key=abs))
ax.set_xlim(xmin=-xabs_max, xmax=xabs_max)
ax.set_box_aspect(1)
ax.set_title('Telescope view of the object')


plt.show()      
    
            

    
