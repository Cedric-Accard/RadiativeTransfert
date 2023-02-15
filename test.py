# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 11:54:00 2023

@author: Cédric
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 14:32:09 2023

@author: Cédric
"""

import numpy as np

pos=np.zeros((1,3))
pos2=np.zeros((1,3))
midpos=np.zeros((1,3))

dire=np.zeros((1,3))

r,r1,r2=0,0,0

R=200.  #Radius of the sphere
norm = 0
count=0
 
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
    for i in range (0,3):
        pos[0,i]=pos[0,i]+dire[0,i]
    ts=ts+(n*sigma + n*sigma)*dl/2
    norm=np.sqrt(pos[0,0]**2+pos[0,1]**2+pos[0,2]**2)

ts2=ts-(n*sigma + n*sigma)*dl/2
count = 1
for i in range (0,3):
    pos2[0,i]=pos[0,i]-dire[0,i]
while (ts-tau>0.001*tau) and (tau-ts2>0.001*tau) :
    for i in range (0,3):
        midpos[0,i]=pos2[0,i]+(1/(2*count))*dire[0,i]        
    count += 1
    x=(ts+ts2)/2
    if x<=tau:
        ts2=x
        for i in range (0,3):
            pos2[0,i]=midpos[0,i]
    else :
        ts=x
        for i in range (0,3):
            pos[0,i]=midpos[0,i]
            
            

    
