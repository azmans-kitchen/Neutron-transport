# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 06:46:45 2020

@author: azmanrafee
"""

import numpy as np

X=[2]
NX=[100]
x0=1
xn=1
azim_div=64
energy_group=6
cores=4

def sigma_scatter(x):
    
    a=np.array([[0.024,0    ,0,0,0    ,0    ],
                [0.171,0.6  ,0,0,0    ,0    ],
                [0.033,0.275,2,0,0    ,0    ],
                [0    ,0    ,0,2,0.275,0.033],
                [0    ,0    ,0,0,0.6  ,0.171],
                [0    ,0    ,0,0,0    ,0.024]])
    return a

def sigma_fis(x):
    a=np.array([0.006,0.06,0.9,0.9,0.06,0.006])
    return a

def sigma_trans(x):
    a=np.array([0.24,0.975,3.1,3.1,0.975,0.24])
    return a

def chi(x):
    a=np.array([0.48,0.02,0,0,0.02,0.48])
    return a

def nu(x):
    a=np.array([3,2.5,2,2,2.5,3])
    return a