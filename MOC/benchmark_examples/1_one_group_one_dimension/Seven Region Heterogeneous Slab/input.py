# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 08:27:16 2020

@author: azmanrafee
"""


""" input description for the case "Abs in Pos. 5" """

f_l=1/0.415                                                     #length of a fuel section
r_l=1/0.371                                                     #length of a reflector section

X=[r_l,f_l,r_l,f_l,r_l,f_l,r_l]                                 #lengths of each section in sequence
NX=[10,10,10,10,10,10,10]                                       #number of meshes for each of the above sections
x0=0                                                            #boundary conditions on the left hand side
xn=0                                                            #boundary conditions on the right hand side
azim_div=64                                                     #number of azimuthal divisions from 0 to 2π
nu=1                                                            #average number of neutrons emitted per fission

def geometry(x):
    f_l=1/0.415
    r_l=1/0.371
    one=(x>0) and (x<r_l)
    two=(x>r_l) and (x<(f_l+r_l))
    three=(x>(f_l+r_l)) and (x<(f_l+(2*r_l)))
    four=(x>(f_l+(2*r_l))) and (x<((2*f_l)+(2*r_l)))
    five=(x>((2*f_l)+(2*r_l))) and (x<((2*f_l)+(3*r_l)))
    six=(x>((2*f_l)+(3*r_l))) and (x<((3*f_l)+(3*r_l)))
    seven=(x>((3*f_l)+(3*r_l))) and (x<((3*f_l)+(4*r_l)))
    fuel=two or four or six
    reflector=one or three or seven
    absorber=(not fuel) and (not reflector)
    
    return fuel,reflector,absorber

def sigma_scatter(x):                                           #scattering cross section as a function of position
    fuel,reflector,absorber=geometry(x)
    return (fuel*0.334)+(reflector*0.334)+(absorber*0.037)

def sigma_trans(x):                                             #(transport corrected) total cross section as a function of position
    fuel,reflector,absorber=geometry(x)
    return (fuel*0.415)+(reflector*0.371)+(absorber*0.371)

def sigma_fis(x):                                               #fission cross section as a function of position
    fuel,reflector,absorber=geometry(x)
    return (fuel*0.178)+(reflector*0)+(absorber*0)