# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 03:37:37 2019

@author: azmanrafee
"""

def moc1Dmg():
    
    import math, numpy as np
    from input import X,NX,x0,xn,azim_div,energy_group,sigma_trans,sigma_scatter,sigma_fis,chi,nu
    
    delx=np.array([])
    N=sum(NX)
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
    phi_j=np.zeros((len(azim),energy_group,j))
    phi_j_boun=np.zeros((2,len(azim),j))
    phi_j_last=np.random.randint(1,10,np.shape(phi_j_boun))/10
    
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
    q_last=q.copy()
    Q=q.copy()
    F_source=np.ones((energy_group,N))
    f_source=F_source.copy()
    S_source=np.ones((energy_group,N))
    for i in range(N):
        S_source[:,i]=np.dot(np.tril(Sigma_scatter[i],k=-1)+np.triu(Sigma_scatter[i],k=1),phi[:,i])
        F_source[:,i]=np.sum(Nu[:,i]*Sigma_fis[:,i]*phi[:,i])*Chi[:,i]
    q_ext=F_source+S_source
    q_ext_last=q_ext.copy()

    iteration=0
    while True:
        iteration+=1
        print(iteration)

        

        
       
        for group in range(energy_group):  
                                
            while True:
                for i in range(N):
                    q[group,i]=q_ext[group,i]+(Sigma_scatter[i,group,group]*phi[group,i])
                    
                Q[group,:]=q[group,:]/(4*math.pi)
                phi[group,:]=q[group,:]/Sigma_trans[group,:]
                
                for angle in range(len(azim)):
                    
                    position=0
                    
                    reverse=False
                    
                    while True:
                        
                        for jj in range(j):
                    
                            s=delx[position]/sintheta_j[jj]/math.cos(azim[angle])
                            
                            del_phi=(phi_j[angle,group,jj]-(Q[group,position]/Sigma_trans[group,position]))*(1-math.exp(-s*Sigma_trans[group,position]))
                            
                            phi[group,position]+=omega_j[jj]*omega_azim*del_phi/(s*Sigma_trans[group,position])
                            
                            phi_j[angle,group,jj]-=del_phi
                        
                        position+=(-1)**reverse
                        if position==N:
                            reverse=True
                            position-=1
                            phi_j_boun[0,angle]=phi_j[angle,group]
                            phi_j[angle,group]*=xn
                        elif position==-1:
                            phi_j_boun[1,angle]=phi_j[angle,group]
                            phi_j[angle,group]*=x0
                            break

                res=np.array([np.max(abs(phi_j_boun-phi_j_last)/phi_j_boun),np.max(abs(phi[group,:]-phi_last[group,:])/phi[group,:])])
                if np.max(res)<1e-8:
                    break
                else:
                    phi_j_last=phi_j_boun.copy()
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
               

    return keff,phi

import time
initial_time=time.time()
keff,phi=moc1Dmg()
print(time.time()-initial_time)
