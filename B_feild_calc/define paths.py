# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 15:58:40 2017

@author: Michael
"""
import numpy as np

#define paths
h=49.1

p1=np.array([[-20.5,25,0-h],
    [19,25,0-h],
    [19,25,104.5-h],
    [-20.5,25,104.5-h],
    [-20.5,25,0-h]])
I1=1
    
p2=np.array([[-20,-26.7,0-h],
    [19.5,-26.7,0-h],
    [19.5,-26.7,108.5-h],
    [-20,-26.7,108.5-h],
    [-20,-26.7,0-h]])
I2=I1    
p3=np.array([[-26,20,0-h],
    [-26,-19,0-h],
    [-26,-19,108.5-h],
    [-26,20,108.5-h],
    [-26,20,0-h]])
I3=1

p4=np.array([[25.2,18.5,0-h],
    [25.2,-20,0-h],
    [25.2,-20,108.5-h],
    [25.2,18.5,108.5-h],
    [25.2,18.5,0-h]])
I4=I3

paths_all=[p1/100,p2/100,p3/100,p4/100]
paths1=[p1/100,p2/100]
paths2=[p3/100,p4/100]

I_all=[I1,I2,I3,I4]
I1=[I1,I2]
I2=[I3,I4]

prev=[1,1,1]

for paths in [paths1,paths2]:
    for p in paths:
        prev=[1,1,1]
        p=p.tolist()
        for point in p:
            if prev==[1,1,1]:
                prev=point
                print()
                continue
            else:  
                print("linecurrent(\"wcs\", \"I\", {0},{1},{2}, {3},{4},{5}, I1);".format(*(*prev,*point)))
                prev=point
        
        
