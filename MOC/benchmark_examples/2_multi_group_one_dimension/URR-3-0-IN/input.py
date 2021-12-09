# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 03:21:07 2019

@author: azmanrafee
"""


import numpy as np

X=[2]
NX=[100]
x0=1
xn=1
azim_div=64
energy_group=3          #total number of energy meshes
cores=4                 #set the number of cores in your system you intend to use


def sigma_scatter(x):
    
    a=np.array([[0.024,0,0],[0.171,0.6,0],[0.033,0.275,2]])
    return a

def sigma_fis(x):
    a=np.array([0.006,0.06,0.9])
    return a

def sigma_trans(x):
    a=np.array([0.24,0.975,3.1])
    return a

def chi(x):
    a=np.array([0.96,0.04,0])
    return a

def nu(x):
    a=np.array([3,2.5,2])
    return a