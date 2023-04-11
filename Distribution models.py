# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 00:24:53 2023

@author: Cédric
"""
import numpy as np
R=200
#Various dnesity models

## Linear
n=np.zeros(R+int(0.1*R)) 
nc=100                      # Central Density number
for j in range(0,R):
    n[j]=nc*(1-j/R)         # Linear model decreasing with R, whith n(R)=0
##


## Onion Layer
n=np.zeros(R+int(0.1*R))
for i in range (0,int(0.1*R)):
    n[i]=100
for i in range (int(0.1*R)+1, int(0.5*R)) :
    n[i]=50
for i in range(int(0.5*R)+1, R):
    n[i]=10
##


## Gaussian

##
#    