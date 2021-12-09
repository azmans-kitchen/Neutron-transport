# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 02:27:15 2019

@author: azmanrafee
"""

def moc1D1g():
    
    import math, numpy as np
    from input import X,NX,x0,xn,azim_div, nu, sigma_scatter,sigma_trans,sigma_fis
    
    delx=np.array([])
    N=sum(NX)
    for i in range(len(X)):
        delx=np.append(delx,np.ones(NX[i])*X[i]/NX[i])
    

    sintheta_j=np.array([0.166648,0.537707,0.932954])
    omega_j=np.array([0.046233,0.283619,0.670148])*2
    j=3
    azim=np.arange(np.pi/azim_div,np.pi/2,2*np.pi/azim_div)   
    omega_azim=4*math.pi/azim_div

    phi=np.ones(N)
    phi_last=phi.copy()
    phi_j=np.zeros((len(azim),j))
    phi_j_boun=np.zeros((2,len(azim),j))
    phi_j_last=np.random.randint(1,10,np.shape(phi_j_boun))/10


    keff_last=1
    keff=1
    Sigma_trans=np.ones(N)
    Sigma_scatter=np.ones(N)
    Sigma_fis=np.ones(N)
    for i in range (N):
        x=np.sum(delx[0:i])+(0.5*delx[i])
        Sigma_trans[i]=sigma_trans(x)
        Sigma_scatter[i]=sigma_scatter(x)
        Sigma_fis[i]=sigma_fis(x)
    F_source=nu*Sigma_fis*phi
    q=F_source+(Sigma_scatter*phi)
    
    Q=q/(4*math.pi)
    iteration=0

    while True:
        iteration+=1
        print(iteration)
        while True:

            phi=q/Sigma_trans

            position=0
            
            reverse=False
            
            while True:
                
                for angle in range(len(azim)):
                    
                    for jj in range(j):
                    
                        s=delx[position]/sintheta_j[jj]/math.cos(azim[angle])

                        del_phi=(phi_j[angle,jj]-(Q[position]/Sigma_trans[position]))*(1-math.exp(-s*Sigma_trans[position]))
                        
                        phi[position]+=omega_j[jj]*omega_azim*del_phi/(s*Sigma_trans[position])
                        
                        phi_j[angle,jj]-=del_phi
                    
                position+=(-1)**reverse
                if position==N:
                    reverse=True
                    position-=1
                        
                    phi_j_boun[0]=phi_j
                    phi_j*=xn

                    
                elif position==-1:
                    phi_j_boun[1]=phi_j
                    phi_j*=x0
                    break
            res=np.array([np.max(abs(phi_j_boun-phi_j_last)/phi_j_boun),np.max(abs(phi-phi_last)/phi)])
            if np.max(res)<1e-8:
                break
            else:
                phi_j_last=phi_j_boun.copy()
                phi_last=phi.copy()

        f_source=nu*Sigma_fis*phi
    
        keff=np.sum(f_source*delx)*keff_last/np.sum(F_source*delx)
        print(keff)
        
        if abs(keff-keff_last)<1e-8:
            break
        else:
            keff_last=keff
            phi_last=phi.copy()
            F_source=f_source.copy()
            q=(F_source/keff+(Sigma_scatter*phi))
            Q=q/(4*math.pi)

        

    return keff,phi
            
keff,phi=moc1D1g()
