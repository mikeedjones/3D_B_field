# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 11:59:15 2018

@author: Michael
"""

import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import bio_sav_law as bio
import pdb


def loop(R,start,end,N,z):
    
    coords=[]
    
    for i in range(0,N):
        theta=start+i*((end-start)/N)
        coords.append([theta,R,z])
    
    return coords

def solenoid(R,z0,dz,dr,nz,nr,N,l_join):
    helix=[]
    zj=list(range(0,nz))
    ri=list(range(0,nr))
    for i in ri:
        for j in zj:
            l=loop(R+dr*i,0,2*np.pi-l_join/(2*np.pi*(R+dr*i)),N,z0+dz*j)
            helix.extend(l)
        zj.reverse()
    
#    print(helix)
#    pdb.set_trace()
    
    return np.array(helix)

def polar_2_cart(polar):
    cart=[]
    for p in polar:
        x=p[1]*np.cos(p[0])
        y=p[1]*np.sin(p[0])
        z=p[2]
        cart.append([x,y,z])
        
    return np.array(cart)

#val=[]
#
#B=0
#Bslicex=0
#flag=True
#B_con=[]
#
#for N in range(10,5000,500):
#    h=solenoid(2,0,0,0.0011,1,1,N,0)
#    hc=polar_2_cart(h)
#    ht=hc.transpose()
#    
#    #fig = plt.figure()
#    #ax = fig.gca(projection='3d')
#    #ax.plot(ht[0,:],ht[1,:],ht[2,:])
#    #plt.show()
#
#    
#    B=bio.trapping_field(hc,n=N)
#    z=B[4]
#    Bslicex=np.einsum('...l,...l',B[0],B[0])[:,0,:]
#
#    B_con.append(B[0])










