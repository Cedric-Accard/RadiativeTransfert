# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 22:38:41 2023

@author: Cédric
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 1

@author: Cédric & Jhanene
"""

import numpy as np
import matplotlib.pyplot as plt
from skspatial.objects import Sphere
import time

plt.close('all')

photons=3

def Source(R, origin, densitymodel) :
    start_time = time.time()

    print("Doing source of origin :",origin, "and radius R=",R)
    
    if densitymodel=='constant' : 
        ## Constant
        n=np.zeros(R+int(0.1*R)) 
        for i in range(0,R):
            n[i]=1000
        ##
    elif densitymodel=='linear' :
        ## Linear distribution ##
        n=np.zeros(R+int(0.1*R)) 
        nc=1000                      # Central Density number
        for j in range(0,R):
            n[j]=nc*(1-j/R)         # Linear model decreasing with R, whith n(R)=0
        ## ###### ############ ##
    elif densitymodel=='layer' :
        ## Layer Model
        n=np.zeros(R+int(0.1*R)) 
        for i in range(0,int(0.2*R)) :
            n[i]=1000
        for i in range(int(0.2*R)+1,int(0.5*R)):
            n[i]=918
        for i in range(int(0.5*R)+1,int(0.8*R)):
            n[i]=700
        for i in range(int(0.8*R)+1,R): 
            n[i]=1.293
    
    ## Arrays for keeping track of all partciles positions
    x_trajectories=np.zeros(shape=(photons,int(2*R)))
    y_trajectories=np.zeros(shape=(photons,int(2*R)))
    z_trajectories=np.zeros(shape=(photons,int(2*R)))
    
    #x_trajectories = np.where(x_trajectories==0, np.nan, x_trajectories)
    x_trajectories[x_trajectories==0]=np.nan
    y_trajectories[y_trajectories==0]=np.nan
    z_trajectories[z_trajectories==0]=np.nan

    ## Arrays for telescope projection 
    x_image, z_image, absorbed = np.zeros(photons), np.zeros(photons), np.zeros(photons) 
    
    
    for j in range(photons) :   # Propagating loop to be done for each photon
        start_time2 = time.time()
    
        print('Doing photon : ', j+1)
        ## (Re-)Initialization of arrays and parameters used for each photon separately ##
        pos=np.zeros((1,3))
        for k in range(0,3):
            pos[0,k]=origin[k]
        pos2=np.zeros((1,3))
        rightpos=np.zeros((1,3))
        dire=np.zeros((1,3))
    
        final_pos=np.zeros((1,3))
        #final_dir=np.zeros((1,3))
        
        norm=0
        count = 0
            # Initialization of arrays to store future values of positions and add starting point with the source origin
        x2, y2, z2, norm_arr = np.array([]), np.array([]), np.array([]), np.array([])
        x2, y2, z2 =np.append(x2,origin[0]), np.append(x2,origin[1]), np.append(x2,origin[2])  
        
        ## ################### ## ###### ### ########## #### ### #### ##### ########### ##
        lambd = float(np.random.uniform(low=400e-9, high=800e-9, size=(1)))
        energy=6.626e-34 * 3e8 / lambd
        e0=6.626e-34 * 3e8 / 400e-9
        ratio = energy/e0
        
        while (norm<=R) :       # While the photon is inside the oject : propagation with successive scattering
        
            #Random number generators for each scattering of each particle
            r, r1, r2 = float(np.random.rand(1)), float(np.random.rand(1)), float(np.random.rand(1))
            abso = float(np.random.rand(1))
            abso*=(1-ratio)
            sigma=0.3           # Cross-Section
    
            tau=-np.log(1-r)    # Optical Depth
            ts=0    
        
            # Definition of the direction vector    
            phi=2*np.pi*r1
            theta=np.pi*r2
            dire[0,0]=np.sin(theta)*np.cos(phi)
            dire[0,1]=np.sin(theta)*np.sin(phi)
            dire[0,2]=np.cos(theta)
           
            while ts<tau :      # Until the ts reach the optical depth propagate forward
                
                if norm<=R :    # Going forward by adding direction vector to the postion while inside the object
                    for i in range (0,3):
                        pos[0,i]=pos[0,i]+dire[0,i]
                    norm=np.sqrt((pos[0,0]-origin[0])**2+(pos[0,1]-origin[1])**2+(pos[0,2]-origin[2])**2)
                    norm_arr=np.append(norm_arr,norm)
                    N=int(norm)
                    dl=tau/(n[0]*sigma)/100 #Integration path
                    ts=ts+ (n[N+1]*sigma + n[N]*sigma)*dl/2                    
                else :
                    break       # Stopping condition if the photon is going out of the object
                count += 1
            
            if norm<=R :        # Correction of the overshoot if inside the object
                ts2=ts-(n[N+1]*sigma + n[N]*sigma)*dl/2   
                for i in range (0,3):
                    pos2[0,i]=pos[0,i]-dire[0,i]
                    
                diffts=ts-ts2
                ratio=(tau-ts2)/diffts
                    
                                # Correction to the position to get forward as close as possible to the limit of the object
                for i in range (0,3): 
                    rightpos[0,i]=pos2[0,i]+ratio*dire[0,i]              
                        
                x2=np.append(x2,rightpos[0,0])  # Adding the position to a list to keep track of them
                y2=np.append(y2,rightpos[0,1])  #
                z2=np.append(z2,rightpos[0,2])  #
            
                if abso>1 :
                    absorbed[j]=1
                    break
            
            else :              # Correction of the overshoot if outside the object
                norm2=norm_arr[count-2]
                for i in range (0,3):
                    pos2[0,i]=pos[0,i]-dire[0,i]    
                
                diff_norm=norm-norm2            # Calculation to get the ratio of dire vector to add to the position to correct it
                ratio_norm=(R-norm2)/diff_norm
            
                for i in range (0,3):
                    rightpos[0,i]=pos2[0,i]+ratio_norm*dire[0,i]    # Correction to the position to get backward as close as possible to the limit of the object
                
                x2=np.append(x2,rightpos[0,0])  # Adding the position to a list to keep track of them
                y2=np.append(y2,rightpos[0,1])  # 
                z2=np.append(z2,rightpos[0,2])  #
                                                
            #final_norm=np.sqrt((rightpos[0,0]-origin[0])**2+(rightpos[0,1]-origin[1])**2+(rightpos[0,2]-origin[2])**2)
                                # Last norm 
        
        final_pos=rightpos      # Last position of the particle at the limit of the object
        #final_dir=dire          # Last direction of the particle at the limit of the object
        
        #theta=np.arccos((final_pos[0,2]-origin[2])/final_norm) # Definition of the theta angle
        #phi=np.arctan2((final_pos[0,1]-origin[1]), (final_pos[0,0]-origin[2]))        
    
        #x_image[j]=final_norm*np.cos(phi)*np.sin(theta) # Definition of the horizontal offset
        x_image[j]=(final_pos[0,0])
        z_image[j]=(final_pos[0,2])           # Z projection of the particle postion 
               
        for k in range(len(x2)):
            x_trajectories[j,k]=x2[k]
        for k in range(len(y2)):
            y_trajectories[j,k]=y2[k]
        for k in range(len(z2)):
            z_trajectories[j,k]=z2[k]
            
        # Remove useless points
        x_image = np.where(absorbed==1, np.nan, x_image)
        z_image = np.where(absorbed==1, np.nan, z_image)
        print("Travel time = ",'{:.2f}'.format((time.time() - start_time2)), "seconds ")

    print("Done")
    print("Total time : ",'{:.2f}'.format((time.time() - start_time)), "seconds" )

    return  origin, x_trajectories, y_trajectories, z_trajectories, x_image, z_image, R

o1, x1, y1, z1, xi1, zi1, R1 = Source(int(6e3), [0,0,0], 'layer')


colors1 = plt.cm.autumn(np.linspace(0, 1, photons))
colors2 = plt.cm.winter(np.linspace(0, 1, photons))
colors3 = plt.cm.summer(np.linspace(0, 1, photons))
# RESULT DISPLAY

fig = plt.figure(figsize=plt.figaspect(2.))


ax = fig.add_subplot(1,1,1 , projection='3d')
ax.grid()
ax.scatter(o1[0], o1[1], o1[2], s=0, color='r', alpha=0.5)
for l in range(0,photons) : 
    ax.plot(x1[l,],y1[l,],z1[l,], color=colors2[l])
sphere1 = Sphere([o1[0], o1[1], o1[2]], R1)
sphere1.plot_3d(ax, alpha=0.1, color='b')

ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=20)
ax.set_zlabel('z', labelpad=20)
ax.set_title('Layer Model')


# =============================================================================
# # Telescope where 0,0 on the graph is the position of the origin point
# ax = fig.add_subplot(1, 2, 2)
# ax.scatter(xi1,zi1,  marker=".", color=colors1)
# ax.scatter(xi2,zi2,  marker=".", color=colors2)
# ax.scatter(xi3,zi3,  marker=".", color=colors3)
# ax.set_title('Telescope projection of sources in xOz plan')
# 
# 
# x_left, x_right = ax.get_xlim()
# y_low, y_high = ax.get_ylim()
# ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))
# =============================================================================

plt.show()      
##    