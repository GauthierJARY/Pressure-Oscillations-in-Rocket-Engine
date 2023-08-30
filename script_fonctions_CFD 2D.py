import numpy as np
import math
import matplotlib.pyplot as plt




def pre_equilibre(u,v,u_star,v_star,dx,dy,dt,nu):

    d2udx2 = (u[0:-2,:]-2*u[:,:]+u[1:,:])/dx**2
    d2udy2 = (u[:,0:-2]-2*u[:,:]+u[:,1:])/dy**2
    ududx = u[:,:]*((u[1:,:]-u[0:-2,:])/(2*dx))
    vdudx = 1/4*(v[0:-2,:]+v[:,:]+v[0:-2,1:]+v[:,1:])*((u[:,1:]-u[:,0:-2])/(2*dy))

    u_star = u + dt * ( nu * (d2udx2 + d2udy2) - (ududx + vdudx) )


    d2vdx2 = (v[0:-2,:]-2*v[:,:]+v[1:,:])/dx**2
    d2vdy2 = (v[:,0:-2]-2*v[:,:]+v[:,1:])/dy**2
    vdvdy = v[:,:]*((v[1:,:]-v[0:-2,:])/(2*dy))
    udvdx = 1/4*(u[0:-2,:]+u[:,:]+u[0:-2,1:]+u[:,1:])*((v[:,1:]-v[:,0:-2])/(2*dx))

    v_star = v + dt * ( nu * (d2vdx2 + d2vdy2) - (udvdx + vdvdy) )




def pression_Poisson(u,v,nx,ny,dx,dy,jmin,jmax,imax,imin):
    # creation of laplacian for poisson equation
    L=np.zeros( nx*ny , nx*ny )
    for j in range(ny):
        for i in range(nx):
            L[i+(j-1)*nx, i+(j-1)*nx] = 2/dx**2+2/dy**2;
            for ii in np.arange(i-1,i+1,2):
                if ( ii >0) & (ii <=nx ): #interior point
                    L[i+(j-1)*nx , ii +(j-1)*nx ]=-1/dx**2
                else: #neumans conditions boundary
                    L[i+(j-1)*nx , i+(j-1)*nx] = L[i+(j-1)*nx , i+(j-1)*nx]-1/dx**2;

            for jj in np.arange(j-1,j+1,2):
                if ( jj >0) &  (jj <=ny ): # interior point
                    L[i+(j-1)*nx,i+(jj-1)*nx]=-1/dy**2;
                else: #Neumans condition
                    L[i+(j-1)*nx,i+(j-1)*nx]= L[i+(j-1)*nx,i+(j-1)*nx]-1/dy**2;
    # pressure in first cells
    L[1,:] = 0
    L[1,1] = 1
    # pressure vector
    n=0
    for j in np.arange(jmin,jmax,1):
        for i in np.arange(imin,imax,1):
            n = n+1
            R[n] = -rho/dt * (   (u_star[i+1,j]-u_star[i,j])/dx + (v_star[i,j+1]-v_star[i,j] )/dy   ) ;
    # pressure matrix
    n=0
    p = np.zeros[imin,imax]
    for j in np.arange(jmin,jmax,1):
        for i in np.arange(imin,imax,1):
            n = n+1
            p[i,j]=pv[n]




def boundary(u,v):
    u[bordures]=0
    v[bordures]=0





def equilibre(u,v,u_star,v_star,p):
    u = u_star - dt*1/rho[:,:]*(p[:,:]-p[0:-1,:])/dx
    v = v_star - dt*1/rho[:,:]*(p[:,:]-p[:,0:-1])/dy


def save(record,t,u,v,p):
    record[0].append(t)
    record[1].append(u)
    record[2].append(v)
    record[3].append(p)