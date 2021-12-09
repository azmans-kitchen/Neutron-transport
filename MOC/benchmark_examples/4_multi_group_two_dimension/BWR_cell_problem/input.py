# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 14:28:40 2020

@author: azmanrafee
"""

import numpy as np

X=8.9
Y=8.9
X0=1
XN=1
Y0=1
YN=1
groups=2
Nx=178
Ny=178
azim_div=40
d=0.05
cores=4

def geometry(x,y):
    fuel=(x>1.5 and x<7.9 and y>1 and y<7.4)
    moderator=not fuel
    return fuel,moderator

def sigma_scatter(x,y):
    
    fuel,moderator=geometry(x,y)
    fuel_cs=np.array([[1.78e-1,1.089e-3],[1.002e-2,5.255e-1]])
    moderator_cs=np.array([[1.995e-1,1.558e-3],[2.188e-2,8.783e-1]])
    return (fuel*fuel_cs)+(moderator*moderator_cs)
    
def sigma_fis(x,y):

    fuel,moderator=geometry(x,y)
    fuel_cs=np.array([6.203e-3,1.101e-1])
    moderator_cs=np.array([0,0])
    return (fuel*fuel_cs)+(moderator*moderator_cs)

def sigma_trans(x,y):

    fuel,moderator=geometry(x,y)
    fuel_cs=np.array([1.96647e-1,5.96159e-1])
    moderator_cs=np.array([2.22064e-1,8.87874e-1])
    return (fuel*fuel_cs)+(moderator*moderator_cs)

def chi(x,y):
    fuel,moderator=geometry(x,y)
    fuel_chi=np.array([1,0])
    moderator_chi=np.array([0,0])
    return (fuel*fuel_chi)+(moderator*moderator_chi)