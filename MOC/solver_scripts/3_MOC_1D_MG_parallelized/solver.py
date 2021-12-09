# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 11:09:46 2019

@author: azmanrafee
"""

def sweepbyangle(energy_group,X0, XN, omega_j, Q_mp, Sigma_trans, omega_azim, NN, delx, sintheta_j,j, group_mp, angles, phi_j_mp, a, b, phi_mpp, on, active):

    import math, numpy as np
    phi_j=np.zeros((energy_group,len(angles),j),dtype='d')
    phi_j_boun=np.zeros((2,len(angles),j),dtype='d')
    phi=np.zeros(NN,dtype='d')
    Q=np.zeros(np.shape(Sigma_trans),dtype='d')
    
    while on.value:
        if active.value:
            phi*=0
            Q[:]=np.asarray(Q_mp).reshape(np.shape(Sigma_trans))
            group=group_mp.value
            for angle in range(len(angles)):
                                                 
                position=0

                              
                reverse=False
                                
                while True:
                    
                    for jj in range(j):
                    
                        s=delx[position]/sintheta_j[jj]/math.cos(angles[angle])
                
                        del_phi=(phi_j[group,angle,jj]-(Q[group,position]/Sigma_trans[group,position]))*(1-math.exp(-s*Sigma_trans[group,position]))
                        
                        phi[position]+=omega_j[jj]*omega_azim*del_phi/(s*Sigma_trans[group,position])                
                        phi_j[group,angle,jj]-=del_phi
        
                                               
                    position+=(-1)**reverse
                    if position==NN:
                        reverse=True
                        position-=1
                        phi_j_boun[0,angle]=phi_j[group,angle]
                        phi_j[group,angle]*=X0
                    elif position==-1:
                        phi_j_boun[1,angle]=phi_j[group,angle]
                        phi_j[group,angle]*=XN
                        break
            phi_mpp[:]=phi
            phi_j_mp[a*2*j:b*2*j]=phi_j_boun.reshape(len(angles)*2*j)
            active.value=False


def moc1Dmg():
    import multiprocessing as mp
    import math, numpy as np
    from ctypes import c_bool
    if __name__=='__main__':
        from input import sigma_scatter,sigma_fis,sigma_trans,chi,nu,X,NX,x0,xn,azim_div,energy_group,cores
        X0=x0
        XN=xn
        delx=np.array([])
        N=sum(NX)
        NN=N
        for i in range(len(X)):
            delx=np.append(delx,np.ones(NX[i])*X[i]/NX[i])
        sintheta_j=np.array([0.166648,0.537707,0.932954])
        omega_j=np.array([0.046233,0.283619,0.670148])*2
        j=3
        azim=np.arange(np.pi/azim_div,np.pi/2,2*np.pi/azim_div)   
        omega_azim=4*math.pi/azim_div
        
        phi=np.ones((energy_group,N))
        phi_last=phi.copy()
        keff_last=1
        Sigma_trans=np.ones((energy_group,N))
        Sigma_scatter=np.ones((N,energy_group,energy_group))
        Sigma_fis=np.ones((energy_group,N))
        Chi=np.ones((energy_group,N))
        Nu=np.ones((energy_group,N))

        for i in range (N):
            x=np.sum(delx[0:i])+(0.5*delx[i])
            Sigma_trans[:,i]=sigma_trans(x)
            Sigma_scatter[i]=sigma_scatter(x)
            Sigma_fis[:,i]=sigma_fis(x)
            Chi[:,i]=chi(x)
            Nu[:,i]=nu(x)
    
        q=np.ones((energy_group,N))
        F_source=np.ones((energy_group,N))
        S_source=np.ones((energy_group,N))
        for i in range(N):
            S_source[:,i]=np.dot(np.tril(Sigma_scatter[i],k=-1)+np.triu(Sigma_scatter[i],k=1),phi[:,i])
            F_source[:,i]=np.sum(Nu[:,i]*Sigma_fis[:,i]*phi[:,i])*Chi[:,i]
        f_source=F_source.copy()
        q_ext=F_source+S_source
        for i in range(N):
            q[:,i]=q_ext[:,i]+(np.diagonal(Sigma_scatter[i])*phi[:,i])
        
        q_last=q.copy()
        Q=q/(4*math.pi)
        
        
        
        iteration=0
            
        group_mp=mp.RawValue('i')
        Q_mp=mp.RawArray('d',np.size(Q))
        Q_mp[:]=Q.reshape(np.size(Q))
        Q=np.asarray(Q_mp).reshape(np.shape(Sigma_trans))
        on=mp.RawValue(c_bool,True)
        active=[]
        angles=[]
        angle_per_core=len(azim)/cores
        range_list=[]
        for i in range(cores+1):
            a_b=round(i*angle_per_core)
            range_list.append(a_b)
        for core in range(cores):
            active.append(mp.RawValue(c_bool,False))
        for core in range(cores):  
            angles.append(azim[range_list[core]:range_list[core+1]])
    
        phi_mp=[]
        for core in range (cores):
            phi_mp.append(mp.RawArray('d',N))
        phi_j_mp=mp.RawArray('d',len(azim)*2*j)
        phi_j_mp[:]=np.zeros(len(phi_j_mp))
        phi_j_last=np.asarray(phi_j_mp[:])
        phi_j=np.zeros((energy_group,len(angles[cores-1]),j),dtype='d')
        phi_j_boun=np.zeros((2,len(angles[cores-1]),j),dtype='d')
        process_list=[]
        for core in range(cores-1):
            t=mp.Process(target=sweepbyangle,args=(energy_group,X0,XN,omega_j,Q_mp,Sigma_trans,omega_azim,NN,delx,sintheta_j,j,group_mp,angles[core],phi_j_mp,range_list[core],range_list[core+1],phi_mp[core],on,active[core]))
            t.start()
            process_list.append(t)

        while True:
            iteration+=1
            print(iteration)
                
            for group in range(energy_group):
                group_mp.value=group
                while True:
                    Q_mp[:]=q.reshape(np.size(q))/(4*math.pi)
                    phi[group,:]=q[group,:]/Sigma_trans[group,:]

                    
                    for core in range(cores-1):
                        active[core].value=True

                    for angle in range(len(angles[cores-1])):

                        position=0
                        reverse=False
                                        
                        while True:
                            for jj in range(j):

                                s=delx[position]/sintheta_j[jj]/math.cos(angles[cores-1][angle])
                        
                                del_phi=(phi_j[group,angle,jj]-(Q[group,position]/Sigma_trans[group,position]))*(1-math.exp(-s*Sigma_trans[group,position]))
                                
                                phi[group,position]+=omega_j[jj]*omega_azim*del_phi/(s*Sigma_trans[group,position])                
                                phi_j[group,angle,jj]-=del_phi
                
                                                       
                            position+=(-1)**reverse

                            if position==NN:
                                reverse=True
                                position-=1
                                phi_j_boun[0,angle]=phi_j[group,angle]
                                phi_j[group,angle]*=X0
                            elif position==-1:
                                phi_j_boun[1,angle]=phi_j[group,angle]
                                phi_j[group,angle]*=XN
                                break
                    phi_j_mp[range_list[cores-1]*2*j:range_list[cores]*2*j]=phi_j_boun.reshape(len(angles[cores-1])*2*j)
                    while np.sum(np.asarray(active)):
                         pass
                        
                    for core in range(cores):
                        phi[group,:]+=phi_mp[core][:]
                        phi_mp[core][:]=np.zeros(N)
                    for i in range(N):
                        q[group,i]=q_ext[group,i]+(Sigma_scatter[i,group,group]*phi[group,i])
    
                    res=np.array([np.max(abs(phi_j_mp[:]-phi_j_last)/phi_j_mp),np.max(abs(phi[group,:]-phi_last[group,:])/phi[group,:])])
                    if np.max(res)<1e-8:
                        break
                    else:
                        phi_j_last=np.asarray(phi_j_mp[:])
                        phi_last[group,:]=phi[group,:].copy()

            for i in range(N):
                f_source[:,i]=np.sum(Nu[:,i]*Sigma_fis[:,i]*phi[:,i])*Chi[:,i]
            keff=np.sum(f_source)*keff_last/np.sum(F_source)
            print(keff)
            keff_last=keff
            q-=q_ext        
            F_source=f_source.copy()
            for i in range(N):
                S_source[:,i]=np.dot(np.tril(Sigma_scatter[i],-1)+np.triu(Sigma_scatter[i],1),phi[:,i])
            q_ext=(F_source/keff)+S_source
            q+=q_ext
            convergence=np.max(abs((q-q_last))/q)
            
            if convergence<1e-6:
                break
            else:
                q_last=q.copy()
        on.value=False
        return phi
import time
initial_time=time.time()
phi=moc1Dmg()
print(time.time()-initial_time)
