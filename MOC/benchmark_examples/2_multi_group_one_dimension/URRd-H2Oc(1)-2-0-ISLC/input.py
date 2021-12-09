# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 11:30:34 2020

@author: azmanrafee
"""

import numpy as np

X=[0.341011,0.751023]
NX=[100]*2
x0=1
xn=1
azim_div=64
energy_group=2
cores=4

def geometry(x):
    aa=0.341011
    bb=0.751023
    fuel=(x>0) and (x<aa)
    moderator=(x>aa) and (x<(aa+bb))
    return fuel,moderator

def sigma_scatter(x):
    fuel,moderator=geometry(x)
    fuel_cs=np.array([[0,0],[0.0342008,2.0688]])
    moderator_cs=np.array([[1.226381244,0],[0.104639534,4.3547]])
    return (moderator * moderator_cs)+(fuel * fuel_cs )

def sigma_fis(x):
    fuel,moderator=geometry(x)
    fuel_cs=np.array([0.61475,0.045704])
    moderator_cs=np.array([0,0])
    return (moderator * moderator_cs)+(fuel * fuel_cs )

def sigma_trans(x):
    fuel,moderator=geometry(x)
    fuel_cs=np.array([0.650917,2.138])
    moderator_cs=np.array([1.331518007,4.3735])
    return (moderator * moderator_cs)+(fuel * fuel_cs )

def chi(x):
    fuel,moderator=geometry(x)
    fuel_chi=np.array([1,0])
    moderator_chi=np.array([0,0])
    return (moderator * moderator_chi)+(fuel * fuel_chi )

def nu(x):
    fuel,moderator=geometry(x)
    fuel_nu=np.array([1.004,2.5])
    moderator_nu=np.array([1,1])
    return (moderator * moderator_nu)+(fuel * fuel_nu )
    
