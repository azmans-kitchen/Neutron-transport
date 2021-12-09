# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 05:00:11 2019

@author: azmanrafee
"""

import numpy as np

X=10                #length along x axis
Y=10                #length along y axis
X0=0                #albedo on lower horizontal boundary
XN=1                #albedo on upper horizontal boundary
Y0=1                #albedo on left vertical boundary
YN=0                #albedo on right vertical boundary
Nx=100              #number of mesh divisions along x axis
Ny=100              #number of mesh divisions along y axis
azim_div=40         
d=0.05              #spacing between parallel rays 

def geometry(x,y):
    fuel=(x>1 and x<2 and y>1 and y<10) or (x>4 and x<5 and y>1 and y<10) or (x>7 and x<8 and y>1 and y<10)
    moderator=not fuel
    return fuel,moderator

def sigma_scatter(x,y):
    
    fuel,moderator=geometry(x,y)
    fuel_cs=1.35
    moderator_cs=0.93
    return (fuel*fuel_cs)+(moderator*moderator_cs)
    
def sigma_fis(x,y):         # nu * sigma_fission

    fuel,moderator=geometry(x,y)
    fuel_cs=0.24
    moderator_cs=0
    return (fuel*fuel_cs)+(moderator*moderator_cs)

def sigma_trans(x,y):

    fuel,moderator=geometry(x,y)
    fuel_cs=1.5
    moderator_cs=1
    return (fuel*fuel_cs)+(moderator*moderator_cs)