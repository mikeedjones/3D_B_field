# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 15:51:30 2018

@author: Michael
"""

import solenoid_3d as s3d
import solenoid_2d as s2d
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

N=500
n=25

def line(roe,I):
    return np.nan_to_num(sc.mu_0*I/(2*sc.pi*np.abs(roe)))

r2d,rho,z=s2d.grid_2d([0,2],[-1,1],n)

#B2d=s2d.print_fields_2d(*s2d.solenoid_2d_map(2,0,0.0,0.0,1,1,r2d,1));

B2d=s2d.print_fields_2d([s2d.solenoid_2d_map(2,0,0,0,1,1,r2d,1)[0][2]],r2d,[[2,2,1]])

#for N in range(100,101,50):
h=s3d.solenoid(2,0,0,0,1,1,N,0)
hc=s3d.polar_2_cart(h)
ht=hc.transpose()
    
#hc=np.array([[-10,0,0],[10,0,0]])

r3d=s3d.bio.r_grid([0,2],[0,2],[-1,1],n)

B3d=s3d.bio.trapping_field(hc,r3d,n=N)

Bslicex, Bslicey, Bslicez=s3d.bio.print_fields(*B3d,[222,222,222]);

B2d_mag=B2d[0]

diffy=B2d_mag-Bslicex
divy=B2d_mag/Bslicex

#std.append(np.std(div))

t=s2d.print_fields_2d([diffy],r2d,[[2,2,3]])
#plt.colorbar()
t=s2d.print_fields_2d([divy],r2d,[[2,2,4]])
#plt.colorbar()

plt.tight_layout()
