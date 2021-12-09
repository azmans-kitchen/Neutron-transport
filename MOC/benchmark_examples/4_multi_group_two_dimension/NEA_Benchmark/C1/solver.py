# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 00:27:22 2020

@author: azmanrafee
"""

def sweepbyangle(j,on,active,group_mp,exponentials,omega_azim,delr,phi_j_init,Sigma_trans,Ny,Nx,Q_mp,area,phi_mp,omega_j,sintheta_j,phi_j_boun_mp,groups,b_points,b_pointer,X0,XN,Y0,YN):

    import numpy as np
    Q=np.asarray(Q_mp).reshape((Ny,Nx))
    phi=np.asarray(phi_mp).reshape((Ny,Nx))
    phi_j_boun=np.asarray(phi_j_boun_mp).reshape((groups,j,2*b_points))
    b_pointer*=2

    while on.value:

        if active.value:

            phi*=0
            group=group_mp.value
            s_point=0
            b_point=b_pointer

            for angle in range(len(omega_azim)):
    
                rays=exponentials[group][angle][0]
                area_mat=exponentials[group][angle][1]
                omega=omega_azim[angle]
                dr=delr[angle]
                
                for ray in rays:
    
                    xpos,ypos=ray[0]
                    tracker=ray[1]
                    exps=ray[2]
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
                            phi[Ny-ypos,xpos-1]+=(del_phi/(area_mat[orient,Ny-ypos,xpos-1]*s_t))*omega*dr*omega_j[jj]*sintheta_j[jj]
    
    
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

            active.value=False


def moc2Dmg():

    #importing necessary libraries
    import multiprocessing as mp
    from ctypes import c_bool
    import numpy as np
    import math

    if __name__=='__main__':

        #importing input data
        from input import X,Y,X0,XN,Y0,YN,geometry,sigma_scatter,sigma_fis,chi,sigma_trans,groups,cores,Nx,Ny,azim_div,d

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

    
        #initiating flux, cross-sections, eigenvalue, and neutron sources 
    
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
        print("initiating ray tracing")
    
        sequence=[]
        from raytracer2 import traceseq2
    
        for i in range(len(azim)):
            sequence .append(traceseq2(X,Y,Nx,Ny,azim[i],delr[i]))
            
        print("ray tracing completed")    
    
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
    

        #preparing data for parallelization
        group_mp=mp.RawValue('i')
        Q_mp=mp.RawArray('d',Ny*Nx)
        Q=np.asarray(Q_mp).reshape((Ny,Nx))
        on=mp.RawValue(c_bool,True)
        active=[]
        for core in range(cores-1):
            active.append(mp.RawValue(c_bool,False))
        angle_per_core=len(azim)/cores
        range_list=[]
        for i in range(cores+1):
            a_b=round(i*angle_per_core)
            range_list.append(a_b)
        phi_mp=[]
        for core in range(cores):
            phi_mp.append(mp.RawArray('d',Nx*Ny))

        expon_core=[]
        for core in range(cores):
            expon_core_groups=[]
            for group in range(groups):
                exponentials_group=exponentials[group][range_list[core]:range_list[core+1]]
                expon_core_groups.append(exponentials_group)
            expon_core.append(expon_core_groups)

        delr_list=[]
        omega_azim_list=[]
        for core in range(cores):
            delr_list.append(delr[range_list[core]:range_list[core+1]])
            omega_azim_list.append(omega_azim[range_list[core]:range_list[core+1]])

        #preparing and declaring boundary conditions for storage
    
        phi_j_init_list=[]
        s_points_list=[]
        for core in range(cores):
            s_points=0
            for angle in exponentials[0][range_list[core]:range_list[core+1]]:
                s_points+=len(angle[0])
            s_points_list.append(s_points)
            phi_j_init_list.append(np.zeros((groups,2,j,s_points)))
    
        phi_j_init=phi_j_init_list[cores-1]
    
        #preparing and declaring boundary conditions for convergence
    
        b_points=0
        b_points_list=[]
        for core in range(cores):
            b_points_list.append(b_points)
            for angle in sequence[range_list[core]:range_list[core+1]]:
        
                nx=angle[0]
                ny=angle[1]
        
                b_points+=(nx*(not X0==1))+(nx*(not XN==1))+(ny*(not Y0==1))+(ny*(not YN==1))
        
            if X0==1:
                b_points+=s_points_list[core]
            
        b_pointer=b_points_list[cores-1]*2

        phi_j_last=np.random.randint(1,10,(groups,j,2*b_points))/10
        phi_j_boun_mp=mp.RawArray('d',np.size(phi_j_last))
        phi_j_boun=np.asarray(phi_j_boun_mp).reshape((groups,j,2*b_points))

        process_list=[]
        for core in range(cores-1):
            t=mp.Process(target=sweepbyangle,args=(j,on,active[core],group_mp,expon_core[core],omega_azim_list[core],delr_list[core],phi_j_init_list[core],Sigma_trans,Ny,Nx,Q_mp,area,phi_mp[core],omega_j,sintheta_j,phi_j_boun_mp,groups,b_points,b_points_list[core],X0,XN,Y0,YN))
            t.start()
            process_list.append(t)
        exponentials=expon_core[cores-1]

        #intiating outer loop
        iteration=0
        converged=False
    
    
        while not converged:
    
            iteration+=1
            print("Outer Iteration No: ",iteration)
    
            #initiating inner loop
    
    
            for group in range(groups):
                group_mp.value=group
                exponentials_group=exponentials[group]
                iter_in=0
                in_converged=False
    
                while not in_converged:
        
                    for ypos in range(Ny):
                        for xpos in range(Nx):
                            q[group,ypos,xpos]=q_ext[group,ypos,xpos]+(Sigma_scatter[ypos,xpos,group,group]*phi[group,ypos,xpos])
    
                    Q[:]=q[group]/(4*math.pi)
                    phi[group]=q[group]/Sigma_trans[group]
                    for core in range(cores-1):
                        active[core].value=True
    
                    b_point=b_pointer
                    s_point=0
                    iter_in+=1
        
                    for angle in range(len(omega_azim_list[cores-1])):
        
                        rays=exponentials_group[angle][0]
                        area_mat=exponentials_group[angle][1]
                        omega=omega_azim_list[cores-1][angle]
                        dr=delr_list[core][angle]
                        
                        for ray in rays:
        
                            xpos,ypos=ray[0]
                            tracker=ray[1]
                            exps=ray[2]
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
        
                    while np.sum(np.asarray(active)):
                        pass

                    for core in range(cores-1):
                        phi[group]+=np.asarray(phi_mp[core]).reshape((Ny,Nx))

                    #print(phi_j_boun[group])
                    res=np.array([np.max(abs((phi[group]-phi_last[group])/phi[group])),np.max(abs((phi_j_boun[group]-phi_j_last[group])/phi_j_boun[group]))])
                    in_converged=np.max(res)<1e-8
                    phi_j_last[group]=phi_j_boun[group].copy()
                    phi_last[group]=phi[group].copy()
                    #print(iter_in,res)
    
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
    
        on.value=False
        return phi

import time
t=time.time()
phi=moc2Dmg()
print("runtime:")
print(time.time()-t)