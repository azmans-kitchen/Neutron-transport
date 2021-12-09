# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:18:13 2020

@author: azmanrafee
"""

import numpy as np

X=96
Y=86
X0=0
XN=0
Y0=0
YN=0
Nx=96
Ny=86
azim_div=32
d=0.5

def geometry(x,y):
    mat1=x>18 and x<48 and y>18 and y<43
    mat2=x>48 and x<78 and y>18 and y<43
    mat3=x>48 and x<78 and y>43 and y<68
    mat4=x>18 and x<48 and y>43 and y<68
    mat5=not (mat1 or mat2 or mat3 or mat4)
    return np.array([mat1,mat2,mat3,mat4,mat5])

def sigma_scatter(x,y):
    
    mat=geometry(x,y)
    cs=np.array([0.53,0.20,0.66,0.50,0.89])
    return np.sum(mat*cs)
    
def sigma_fis(x,y):

    mat=geometry(x,y)
    cs=np.array([0.079,0,0.043,0,0])
    return np.sum(mat*cs)

def sigma_trans(x,y):

    mat=geometry(x,y)
    cs=np.array([0.6,0.48,0.7,0.65,0.9])
    return np.sum(mat*cs)
