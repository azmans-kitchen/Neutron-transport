# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 11:29:31 2020

@author: azmanrafee
"""

def moc2Dmg():

    #importing input data
    from input import X,Y,X0,XN,Y0,YN,sigma_scatter,sigma_fis,chi,sigma_trans,groups,Nx,Ny,azim_div,d

    #importing necessary libraries
    import numpy as np, math

    #mesh size
    delx=X/Nx
    dely=Y/Ny
    area=delx*dely

    #yamamoto takeshi polar quarature
    sintheta_j=np.array([0.166648,0.537707,0.932954]).reshape((3,1))
    omega_j=np.array([0.046233,0.283619,0.670148])*2
    j=3

    #azimuthal quadrature and distances between lines

    angles=np.arange(math.pi/azim_div,math.pi/2,2*math.pi/azim_div)
    delr=np.zeros(len(angles))
    azim=delr.copy()


    for i in range(len(angles)):

        nx=math.ceil(X*abs(math.sin(angles[i]))/d)
        ny=math.ceil(Y*abs(math.cos(angles[i]))/d)
        dx=X/nx
        dy=Y/ny
        dr=dx*dy/math.sqrt(dx**2+dy**2)
        azim[i]=math.atan(dy/dx)
        delr[i]=dr
    
    omega_azim=np.zeros(len(azim))
    for i in range(len(azim)):
        if i==0:
            a=0
            b=(azim[i+1]+azim[i])/2
        elif i==len(azim)-1:
            a=(azim[i]+azim[i-1])/2
            b=math.pi/2
        else:
            a=(azim[i]+azim[i-1])/2
            b=(azim[i+1]+azim[i])/2
        omega_azim[i]=b-a

    #declaring flux, cross-sections, eigenvalue, and neutron sources 

    phi=np.ones((groups,Ny,Nx))
    phi_last=phi.copy()
    keff_last=1
    Sigma_trans=np.zeros((groups,Ny,Nx))
    Sigma_scatter=np.zeros((Ny,Nx,groups,groups))
    Sigma_fis=np.zeros((groups,Ny,Nx))
    F_source=np.zeros((groups,Ny,Nx))
    f_source=np.zeros((groups,Ny,Nx))
    S_source=np.zeros((groups,Ny,Nx))
    q_ext=np.zeros((groups,Ny,Nx))
    q=np.zeros((groups,Ny,Nx))
    Chi=np.zeros((groups,Ny,Nx))

    #computing cross-sections and initial source

    for xpos in range(Nx):
        x=(xpos+0.5)*delx
        for ypos in range(Ny):
            y=(ypos+0.5)*dely
            Sigma_trans[:,Ny-1-ypos,xpos]=sigma_trans(x,y)
            Sigma_scatter[Ny-1-ypos,xpos]=sigma_scatter(x,y)
            Sigma_fis[:,Ny-1-ypos,xpos]=sigma_fis(x,y)
            Chi[:,Ny-1-ypos,xpos]=chi(x,y)

    for ypos in range(Ny):
        for xpos in range(Nx):
            F_source[:,ypos,xpos]=np.sum(Sigma_fis[:,ypos,xpos]*phi[:,ypos,xpos])*Chi[:,ypos,xpos]
            S_source[:,ypos,xpos]=np.dot(np.tril(Sigma_scatter[ypos,xpos],k=-1)+np.triu(Sigma_scatter[ypos,xpos],k=1),phi[:,ypos,xpos])
    q_ext=F_source+S_source
    q_ext_last=q_ext.copy()
    
    #obtaining ray tracing data

    sequence=[]
    from raytracer2 import traceseq2

    for i in range(len(azim)):
        sequence .append(traceseq2(X,Y,Nx,Ny,azim[i],delr[i]))


    #precomputing exponentials
    print("precomputing exponentials")
    exponentials=[]
    for group in range(groups):
        exponentials_group=[]
    
        for angle in range(len(azim)):
    
            rays=sequence[angle][2]
            norm_mat=sequence[angle][3]
            exp_rays=[]
    
            for ray in rays:
    
                lengths=ray[1]
                exp_ray=np.zeros((j,len(lengths)))
                pointer=0
                tracks=[]
                xpos,ypos=ray[0]
                xreverse=False
                yreverse=False
                orient=int(xreverse==yreverse)
    
                for pr_length in lengths:
    
                    length=pr_length/sintheta_j
                    expon=1-np.exp(-abs(length)*Sigma_trans[group,Ny-ypos,xpos-1])
                    track=pr_length
                    if track:
                        track/=abs(pr_length)
                    else:
                        track=1.0
    
                    exp_ray[:,pointer]=expon[:,0]
                    pointer+=1
    
                    if track==-1:
                        xpos+=(-1)**xreverse
                        tracks+=[False]
                    if track==1 or track==0:
                        ypos+=(-1)**yreverse
                        tracks+=[True]
    
                    if xpos==Nx+1:
                        xpos-=1
                        xreverse=True
                    if xpos==0:
                        xpos+=1
                        xreverse=False
                    if ypos==Ny+1:
                        ypos-=1
                        yreverse==True
                    if ypos==0:
                        ypos+=1
                        yreverse==False
    
                exp_rays.append([ray[0],tracks,exp_ray])
            exponentials_group.append([exp_rays,norm_mat])
        exponentials.append(exponentials_group)
    print("exponential calculations completed")


    #preparing and declaring boundary conditions for storage

    s_points=0

    for angle in exponentials[0]:
        s_points+=len(angle[0])

    phi_j_init=np.zeros((groups,2,j,s_points))

    #preparing and declaring boundary conditions for convergence

    b_points=0

    for angle in sequence:

        nx=angle[0]
        ny=angle[1]

        b_points+=(nx*(not X0==1))+(nx*(not XN==1))+(ny*(not Y0==1))+(ny*(not YN==1))

    if X0==1:
        b_points+=s_points

    phi_j_last=np.random.randint(1,10,(groups,j,2*b_points))/10
    phi_j_boun=np.zeros(np.shape(phi_j_last))


    #intiating outer loop
    iteration=0
    converged=False


    while not converged:

        iteration+=1
        print("Outer Iteration No: ",iteration)

        #initiating inner loop


        for group in range(groups):
            exponentials_group=exponentials[group]
            in_converged=False

            while not in_converged:
    
                for ypos in range(Ny):
                    for xpos in range(Nx):
                        q[group,ypos,xpos]=q_ext[group,ypos,xpos]+(Sigma_scatter[ypos,xpos,group,group]*phi[group,ypos,xpos])

                Q=q[group]/(4*math.pi)
                phi[group]=q[group]/Sigma_trans[group]

                b_point=0
                s_point=0

    
                for angle in range(len(azim)):
    
                    rays=exponentials_group[angle][0]
                    area_mat=exponentials_group[angle][1]
                    omega=omega_azim[angle]
                    dr=delr[angle]
                    
                    for ray in rays:
    
                        xpos,ypos=ray[0]
                        tracker=ray[1]
                        exps=ray[2]
                        #print(s_point,np.shape(phi_j_init))
                        phi_j=phi_j_init[group,0,:,s_point]
                        pointer=0
                        reverse=False
                        xreverse=False
                        yreverse=False
                        orient=int(xreverse==yreverse)
                        points=np.size(exps)/j
    
    
                        while True:
    
                            s_t=Sigma_trans[group,Ny-ypos,xpos-1]
                            Qq=Q[Ny-ypos,xpos-1]
    
                            for jj in range(j):
    
                                del_phi=(phi_j[jj]-(Qq/s_t))*exps[jj,pointer]
                                phi_j[jj]-=del_phi
                                phi[group,Ny-ypos,xpos-1]+=(del_phi/(area_mat[orient,Ny-ypos,xpos-1]*s_t))*omega*dr*omega_j[jj]*sintheta_j[jj]
    
    
                            if tracker[pointer]:
                                ypos+=(-1)**yreverse
                            else:
                                xpos+=(-1)**xreverse
    
                            if xpos==Nx+1:
                                xpos-=1
                                xreverse=True
                                orient=int(xreverse==yreverse)
                                if not YN==1:
                                    phi_j_boun[group,:,b_point]=(1-YN)*phi_j
                                    b_point+=1
                                    phi_j*=YN
    
                            if ypos==Ny+1:
                                ypos-=1
                                yreverse=True
                                orient=int(xreverse==yreverse)
                                if not XN==1:
                                    phi_j_boun[group,:,b_point]=(1-XN)*phi_j
                                    b_point+=1
                                    phi_j*=XN
    
                            if ypos==0:
                                ypos+=1
                                yreverse=False
                                orient=int(xreverse==yreverse)
                                if not X0==1:
                                    phi_j_boun[group,:,b_point]=(1-X0)*phi_j
                                    b_point+=1
                                    phi_j*=X0
    
                            if xpos==0:
                                xpos+=1
                                xreverse=False
                                orient=int(xreverse==yreverse)
                                if not Y0==1:
                                    phi_j_boun[group,:,b_point]=(1-Y0)*phi_j
                                    b_point+=1
                                    phi_j*=Y0
    
                            pointer+=(-1)**reverse
    
                            if pointer==points:
                                reverse=True
                                if X0==1:
                                    phi_j_boun[group,:,b_point]=phi_j
                                    b_point+=1
                                phi_j_init[group,0,:,s_point]=phi_j
                                phi_j=phi_j_init[group,1,:,s_point]
                                tracker=[True]+tracker
                                pointer-=1
                                xreverse=True
                                orient=int(xreverse==yreverse)
    
                            if pointer==-1:
                                if X0==1:
                                    phi_j_boun[group,:,b_point]=phi_j
                                    b_point+=1
                                phi_j_init[group,1,:,s_point]=phi_j
                                break
    
                        s_point+=1
    
                res=np.array([np.max(abs((phi[group]-phi_last[group])/phi[group])),np.max(abs((phi_j_boun[group]-phi_j_last[group])/phi_j_boun[group]))])
                in_converged=np.max(res)<1e-8
                phi_j_last[group]=phi_j_boun[group].copy()
                phi_last[group]=phi[group].copy()


        for ypos in range(Ny):
            for xpos in range(Nx):
                f_source[:,ypos,xpos]=np.sum(Sigma_fis[:,ypos,xpos]*phi[:,ypos,xpos])*Chi[:,ypos,xpos]
                S_source[:,ypos,xpos]=np.dot(np.tril(Sigma_scatter[ypos,xpos],k=-1)+np.triu(Sigma_scatter[ypos,xpos],k=1),phi[:,ypos,xpos])
        keff=np.sum(f_source)*keff_last/np.sum(F_source)
        q_ext=(f_source/keff)+S_source

        print(keff)
        if abs(keff-keff_last)<1e-5:
            converged=True

        keff_last=keff
        q_ext_last=q_ext.copy()
        F_source=f_source.copy()

    return phi

import time
t=time.time()
phi=moc2Dmg()
print("runtime:")
print(time.time()-t)