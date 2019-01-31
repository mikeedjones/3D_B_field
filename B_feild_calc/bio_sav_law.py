# -*- coding: utf-8 -*-
# defines the Bio savat law 

from tqdm import tqdm
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
import pdb
from mpl_toolkits.mplot3d import Axes3D

n=200

def bs(dlhat,rprime):
#    pdb.set_trace()
    return np.divide(np.cross(dlhat, rprime),(np.sqrt(np.einsum('...l,...l',rprime,rprime))**3)[:,:,:,None]);
    
def integrate_bs(l,r,n):
#    pdb.set_trace()
    Batr=np.zeros((3)) 
    for dl, dlhat, I in tqdm(l,total=n):
        rprime = r-dl
        B_temp=I*bs(dlhat,rprime)
        Batr = Batr+B_temp

#        plt.imshow(np.einsum('...l,...l',rprime,rprime)[:,:,12]**3)
#        plt.show()
    
    return sc.mu_0/(4*sc.pi)*Batr;

#the paths are split into n sections of equal length

def split_paths(paths,n,path_currents):
    l_split=[]#vector of the path
    l_split_hat=[] #position vector of the path
    I=[] #current in that section of path
    
    fig=plt.figure()
    
    for p,i in zip(paths,path_currents):
        split=split_path(p,n,i,fig=fig)
        l_split.extend(split[0])
        l_split_hat.extend(split[1])
        I.extend(split[2])
    
    return zip(l_split, l_split_hat,I)

def split_path(p,n,current, fig=plt.figure()):
    # Prime the pump
    totlen = 0
    for cur,nxt in zip( p, p[1:]):
        totlen = totlen+np.sqrt((cur-nxt).dot(cur-nxt))

#    print(totlen)    
    
    dl_len=totlen/n #divide the total length of the path into equal sections
    l_split=[]
    l_hat_split=[]
    I=[]
#    length=[]
    excess=0
    ex=[]
    for cur,nxt in zip( p, p[1:]):
        section_len=np.sqrt((cur-nxt).dot(cur-nxt))
        dlhat=-(cur-nxt)/section_len
        dl=cur+excess*dlhat
        running = True
        while running:
#            if len(l_split)>1:  
#                length.append(np.sqrt((l_split[-1]-dl).dot(l_split[-1]-dl)))        
            l_split.append(dl)
            l_hat_split.append(dlhat*dl_len)
            I.append(current)
            dl=dl+dl_len*dlhat
            
            if np.sqrt((cur-dl).dot(cur-dl))>section_len:
                excess=np.sqrt((cur-dl).dot(cur-dl))-section_len
                ex.append(excess)
                running = False
#    pdb.set_trace()
#    print(len(l_split)) 
    lt=np.array(l_split).transpose()
    A=[]
    i=0
    for el in lt[0,:]:
        A.append(i)
        i+=1
#    pdb.set_trace()
    C=plt.cm.viridis(np.array(A)/A[-1])
#    plt.plot(ex)
#    plt.show()
##    
##    plt.plot(length)
##    print(np.std(length))
    ax = fig.gca(projection='3d')
    ax.scatter(lt[0,:],lt[1,:],lt[2,:],marker='.', color=C)
    plt.show()
    print(max(np.sqrt(lt[0,:]**2+lt[1,:]**2)))
    
    return [l_split,l_hat_split,I];
        
#defines the gridpoints at which the Bfield is evaluated

def r_grid(x,y,z,n):
    x=np.linspace(x[0], x[1], n[0])
    y=np.linspace(y[0], y[1], n[1])
    z=np.linspace(z[0], z[1], n[2])
    
    r = np.meshgrid(x, y, z)

    return np.transpose(r),x,y,z;

def trapping_field(paths,grid,n=100,I=[1]):
    
    r,x,y,z=grid
    
    l=split_paths(paths,n,I)
#    print(l)

    
    B=integrate_bs(l,r,n)
    
    B=np.nan_to_num(B)
    
    return B,r,x,y,z #-np.mean(B[int(y.size/2),int(x.size/2)],axis=0)

def print_fields(B,r,x,y,z,fig,vmax=99):    
    #einsum is used to sum the components of the B field vectors at each point
    
    Bre=B.reshape(x.size*y.size*z.size,3)
    
    Rre=r.reshape(x.size*y.size*z.size,3)
    
    plt.subplot(fig[0])
    Bslicez=np.sqrt(np.einsum('...l,...l',B,B)[int(z.size/2),:,:])
    
    mi=np.percentile(Bslicez.flatten(),[1,vmax])

    plt.imshow(Bslicez, vmin=mi[0], vmax=mi[1],
#               interpolation='mitchell', 
               cmap='viridis', extent=[min(x),max(x),max(y),min(y)],aspect='auto')
#    plt.colorbar()
    
    plt.subplot(fig[1])
    Bslicey=np.sqrt(np.einsum('...l,...l',B,B)[:,int(x.size/2),:])
    
    mi=np.percentile(Bslicey.flatten(),[1,vmax])
    
    plt.imshow(Bslicey, vmin=mi[0], vmax=mi[1],
#                interpolation='mitchell', 
                cmap='viridis', extent=[min(y),max(y),min(z),max(z)],aspect='auto')
    
#    plt.colorbar()
    
    plt.subplot(fig[2])
    Bslicex=np.sqrt(np.einsum('...l,...l',B,B)[:,:,int(y.size/2)])
    
    mi=np.percentile(Bslicex.flatten(),[1,vmax])
    
    plt.imshow(Bslicex, vmin=mi[0], vmax=mi[1],
#                interpolation='mitchell', 
                cmap='viridis', extent=[min(x),max(x),min(z),max(z)],aspect='auto')
#    plt.colorbar()
    np.savetxt('B_trapping_field.txt', np.hstack((Rre,Bre)), header="x\ty\tz\tBx\tBy\tBz")
    
#    plt.show()
#    plt.clf()
    return Bslicex, Bslicey, Bslicez
    
    
def steering_coils(paths,n=200,I=[1],fig=[231,132,133],vmax=99,fname='B_field_steering.txt'):
    N=[101,101,251]
    r,x,y,z=r_grid([-0.001,0.001],[-0.001,0.001],[0.25,0.5],N)
    
    l=split_paths(paths,n,I)
    plt.show()
#    pdb.set_trace()
    B=integrate_bs(l,r,n*len(paths))
#    B=np.nan_to_num(B)-np.mean(B[int(y.size/2),int(x.size/2),10:300,0:1],axis=0)
    B[:,:,:,0]=B[:,:,:,0]-np.mean(B[int(y.size/2),int(x.size/2),10:300,0],axis=0)
    B[:,:,:,1]=B[:,:,:,1]-np.mean(B[int(y.size/2),int(x.size/2),10:300,1],axis=0)
    
    #einsum is used to sum the components of the B field vectors at each point
    
    Bre=B.reshape(x.size*y.size*z.size,3)
    
    Rre=r.reshape(x.size*y.size*z.size,3)

    plt.subplot(fig[0])
    Bslicez=np.sqrt(np.einsum('...l,...l',B,B)[int(z.size)-1,:,:])
    
    mi=np.percentile(Bslicez.flatten(),[1,vmax])
    
    plt.imshow(Bslicez, vmin=mi[0], vmax=mi[1],
#               interpolation='mitchell', 
               cmap='viridis',aspect='auto', extent=[min(x),max(x),max(y),min(y)])
    plt.colorbar()
    
    plt.subplot(fig[1])
    Bslicey=np.sqrt(np.einsum('...l,...l',B,B)[:,int(x.size/2),:])
    
    mi=np.percentile(Bslicey.flatten(),[1,vmax])
    
    plt.imshow(Bslicey, vmin=mi[0], vmax=mi[1],
#                interpolation='mitchell', 
                cmap='viridis',aspect='auto', extent=[min(y),max(y),max(z),min(z)])
    
    plt.colorbar()
    
    plt.subplot(fig[2])
    Bslicex=np.sqrt(np.einsum('...l,...l',B,B)[:,:,int(y.size/2)])
    
    mi=np.percentile(Bslicex.flatten(),[1,vmax])
    
    plt.imshow(Bslicex, vmin=mi[0], vmax=mi[1],
#                interpolation='mitchell', 
                cmap='viridis',aspect='auto', extent=[min(x),max(x),max(z),min(z)])

    plt.tight_layout()
    
    np.savetxt(fname, np.hstack((Rre,Bre)), header="x\ty\tz\tBx\tBy\tBz")

