# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 11:52:22 2020

@author: azmanrafee
"""

def traceseq2 ( X=4 , Y=4 , Nx=40 , Ny=40 , theta=3.1416/6 , d=1/30 ):

    import math,time
    import numpy as np

    xgrid=X/Nx
    ygrid=Y/Ny
    nx=int(round(X*abs(math.sin(theta))/d))
    ny=int(round(Y*abs(math.cos(theta))/d))

    delx = X /  nx
    dely = Y /  ny
    m=dely/delx
    delr = ( delx * dely ) / math.sqrt( delx ** 2 + dely ** 2)

    starting_points=np.ones(nx,dtype=bool)
    sequence=[]
    norm_mat=np.zeros((2,Ny,Nx))

    for point in range(nx):
        if starting_points[point]:
            x0=(0.5+point)*delx
            y0=0
            xpos=math.floor((x0/xgrid)+1)
            ypos=1

            seqnloc=[(xpos,ypos)]
            lengths=np.zeros([(ny+nx)*Nx*Ny])
            pointer=0
            x=x0
            y=y0
            xreverse=False
            yreverse=False
            orient_=1
    
            while starting_points[point]:

                deltax=(xpos-xreverse)*xgrid-x
                deltay=abs(m*deltax)*((-1)**yreverse)
    
                if abs(deltay)<=abs((ypos-yreverse)*ygrid-y):
                    x+=deltax
                    y+=deltay
                    norm_loc=Ny-ypos,xpos-1
                    xpos+=(-1)**xreverse
                    gridtracker=-1
    
                else:
                    deltay=(ypos-yreverse)*ygrid-y
                    deltax=abs(deltay/m)*((-1)**xreverse)
                    x+=deltax
                    y+=deltay
                    norm_loc=Ny-ypos,xpos-1
                    ypos+=(-1)**yreverse
                    gridtracker=1

                s=math.sqrt(deltax**2+deltay**2)
                if s<1e-16:
                    s+=1e-16
                lengths[pointer]=s*gridtracker
                pointer+=1
                norm_mat[orient_][norm_loc]+=s

                if xpos==Nx+1:
                    xreverse=True
                    xpos-=1

                if ypos==Ny+1:
                    yreverse=True
                    ypos-=1

                if xpos==0:
                    xreverse=False
                    xpos+=1

                if ypos==0:
                    yreverse=False
                    ypos+=1
                    skip_point=int(round((x/delx)-0.5))
                    starting_points[skip_point]=False

                orient_=int(xreverse==yreverse)
            lengths=lengths[0:pointer]
            seqnloc.append(lengths)
            sequence.append(seqnloc)

    norm_mat*=delr

    return [nx,ny,sequence,norm_mat]