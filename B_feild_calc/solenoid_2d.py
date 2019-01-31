# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 15:21:35 2018

@author: Michael
"""

import numpy as np
import scipy.special as ss
import scipy.constants as sc
from tqdm import tqdm
import matplotlib.pyplot as plt
import pdb

def K(rho,z,R,z0):
#    pdb.set_trace()
    k=(4*R*rho)/(np.square(R+rho)+np.square(z-z0))
    k=np.nan_to_num(k)
    k=np.clip(k,0,1)
    return ss.ellipk(k)

def E(rho,z,R,z0):
#    pdb.set_trace()
    k=(4*R*rho)/(np.square(R+rho)+np.square(z-z0))
    k=np.nan_to_num(k)
    k=np.clip(k,0,1)
    return ss.ellipe(k)

def Bz(P,R,z0,I):
#    pdb.set_trace()
#    
    rho=P[0]
    z=P[1]
    
    pfac=sc.mu_0*I/(2*sc.pi)
    
    fac1=1/(((R+rho)**2+(z-z0)**2)**0.5)
    
    fac2=K(rho,z,R,z0)+((R**2-rho**2-(z-z0)**2)/((R-rho)**2+(z-z0)**2))*E(rho,z,R,z0)
    
    return pfac*fac1*fac2


#Bcommon=B0./((r+rho(m))^2+(ztest-z).^2).^0.5; %common terms
#Bz=Bcommon.*(K+E.*(r^2-rho(m)^2-(ztest-z).^2)./((r-rho(m))^2+(ztest-z).^2)); %work out Bz for this loop at each radial pos and add it to the total
#Bztot(m,:)=Bztot(m,:)+Bz;


#Br=Bcommon.*(ztest-z)/rho(m).*(-K+E.*(r^2+rho(m)^2+(ztest-z).^2)./((r-rho(m))^2+(ztest-z).^2));


def Brho(P,R,z0,I):    
    
#    pdb.set_trace()
    
    rho=P[0]
    z=P[1]

    pfac=sc.mu_0*I/(2*sc.pi*rho)
    
    fac1=(z-z0)/(((R+rho)**2+(z-z0)**2)**0.5)
    
    fac2=-K(rho,z,R,z0)+((R**2+rho**2+(z-z0)**2)/((R-rho)**2+(z-z0)**2))*E(rho,z,R,z0)

    return pfac*fac1*fac2

def solenoid_2d_map(R,z0,dz,dr,nz,nr,grid,I):
    Bz_map=0
    Brho_map=0
    wires=[]
    for z in tqdm(np.linspace(z0,z0+dz*(nz-1),nz)):
        for r in np.linspace(R,R+dr*(nr-1),nr):
            Bz_map+=Bz(grid,r,z,I)
            Brho_map+=Brho(grid,r,z,I)
            wires.append([r,z])
    
#    print(Bz_map)
#    print(Brho_map)
    
#    Brho_map[:,0]=np.zeros(Brho_map[:,0].shape)
    
    Bz_map=np.nan_to_num(Bz_map)-np.flip(np.nan_to_num(Bz_map),axis=0)
    Brho_map=np.nan_to_num(Brho_map)+np.flip(np.nan_to_num(Brho_map),axis=0)
    magB_map=np.sqrt(np.square(Bz_map)+np.square(Brho_map))
    magB_map[:,0]=abs(Bz_map[:,0])
#    magB_map-=np.flip(np.nan_to_num(magB_map),axis=0)
    print(max(np.array(wires)[0,:]))
    
    magB_map=np.hstack([np.flip(magB_map,axis=1),magB_map])
    
    return [Bz_map, Brho_map, magB_map], grid

def grid_2d(rho,z,n):
    rho=np.linspace(rho[0], rho[1], num=n[0])
    z=np.linspace(z[0], z[1], num=n[1])
    
    r = np.meshgrid(rho,z)

    return r,rho,z;

def print_fields_2d(maps,r,fig,vmax=99):
    
#    magB=Bz
#    pdb.set_trace()
    
    ext=[r[0][0,0],r[0][-1,-1],r[1][0,0],r[1][-1,-1]]
#    ext=[ext[0],ext[2],ext[1],ext[3]]
    
    for B,f in zip(maps,fig):
        plt.subplot(f)
        mi=np.percentile(B.flatten(),[1,vmax])
        plt.imshow(B, vmin=mi[0], vmax=mi[1],#interpolation='mitchell',
                   cmap='viridis', extent=ext, aspect='auto')
#        plt.colorbar()

#        plt.show()
    plt.tight_layout()

    return maps

#r,rho,z=grid_2d([.015,.02],[-0.05,.05],[200,200])
#
#figs=[plt.subplot(231),plt.subplot(233),plt.subplot(234),plt.subplot(235),plt.subplot(236)]
#
#sol=solenoid_2d_map(0.018,0.0095,0.0011,0.0011,12,12,r,1)
#maps=print_fields_2d(*sol,figs)
#
#figs_grad=[]




