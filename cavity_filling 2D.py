############################
### Fonctions definition ###
############################




import numpy as np
import math
import matplotlib.pyplot as plt




def pre_equilibre(u,v,u_star,v_star,dx,dy,dt,nu):
    d2udx2=u.copy()
    d2udy2=u.copy()
    ududx=u.copy()
    vdudx=u.copy()
    d2vdx2=u.copy()
    d2vdy2=u.copy()
    vdvdy=u.copy()
    udvdx=u.copy()
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            d2udx2[i,j] = (u[i-1,j]-2*u[i,j]+u[i+1,j])/dx**2
            d2udy2[i,j] = (u[i,j-1]-2*u[i,j]+u[i,j+1])/dy**2
            ududx[i,j] = u[i,j]*((u[i+1,j]-u[i-1,j])/(2*dx))
            vdudx[i,j] = 1/4*(v[i-1,j]+v[i,j]+v[i-1,j+1]+v[i,j+1])*((u[i,j+1]-u[i,j-1])/(2*dy))
        
            u_star[i,j] = u[i,j] + dt * ( nu * (d2udx2[i,j] + d2udy2[i,j]) - (ududx[i,j] + vdudx[i,j]) )
        
        
            d2vdx2[i,j] = (v[i-1,j]-2*v[i,j]+v[i+1,j])/dx**2
            d2vdy2[i,j] = (v[i,j-1]-2*v[i,j]+v[i,j+1])/dy**2
            vdvdy[i,j] = v[i,j]*((v[i+1,j]-v[i-1,j])/(2*dy))
            udvdx[i,j] = 1/4*(u[i-1,j]+u[i,j]+u[i-1,j+1]+u[i,j+1])*((v[i,j+1]-v[i,j-1])/(2*dx))
        
            v_star[i,j] = v[i,j] + dt * ( nu * (d2vdx2[i,j] + d2vdy2[i,j]) - (udvdx[i,j] + vdvdy[i,j]) )




def pression_Poisson(u,v,nx,ny,dx,dy,jmin,jmax,imax,imin):
    # creation of laplacian for poisson equation
    ni=nx*ny
    L=np.zeros( (ni,ni) )
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
    R=np.zeros(100)
    for j in np.arange(jmin,jmax,1):
        for i in np.arange(imin,imax,1):
            n = n+1
            R[n] = -rho[i,j]/dt * (   (u_star[i+1,j]-u_star[i,j])/dx + (v_star[i,j+1]-v_star[i,j] )/dy   ) ;
    # pressure matrix
    n=0
    pv = np.zeros((imin,imax))
    for j in np.arange(jmin,jmax,1):
        for i in np.arange(imin,imax,1):
            n = n+1
            p[i,j]=pv[n]




def boundary(u,v,U0,R):
    u[bordures==1]=0
    v[bordures==1]=0
    u[round(nx/2)-R:round(nx/2)+R:,0]=0
    v[round(nx/2)-R:round(nx/2)+R:,0]=-U0




def equilibre(u,v,u_star,v_star,p):
    for i in range(1,nx):
        for j in range(1,ny):
            u[i,j] = u_star[i,j] - dt*1/rho[i,j]*(p[i,j]-p[i-1,j])/dx
            v[i,j] = v_star[i,j] - dt*1/rho[i,j]*(p[i,j]-p[i,j-1])/dy


def save(record,t,u,v,p):
    record[0].append(t)
    record[1].append(u)
    record[2].append(v)
    record[3].append(p)




############################
### Loop ###
############################


import numpy as np
import math
import matplotlib.pyplot as plt

L=1
R=1
U0=1
col = 10
row = 10
nx = 10
ny = 10
dx = L/nx
dy = L/ny
dt = 0.1
viscosity = 0.01
time_vector = 0
time_simulation = 10
jmin=1
imin=1
jmax=nx-1
imax=nx-1
boundary_conditions = None
cavity = np.zeros( (col,row) )
u = np.zeros( (col,row) )
u_star = np.zeros( (col,row) )
v_star = np.zeros( (col,row) )
v = np.zeros( (col,row) )
rho = np.zeros( (col,row) )
p = np.zeros( (col,row) )
nu=viscosity
record=[None,None,None,None]

bordures=np.zeros((nx,ny))
bordures[0,:]=1
bordures[-1,:]=1
bordures[:round(nx/2)-R,0]=1
bordures[round(nx/2)+R:,0]=1
bordures[:,-1]=1


while time_vector < time_simulation:
    time_vector+=dt
    boundary(u,v,U0,R)
    pre_equilibre(u,v,u_star,v_star,dx,dy,dt,nu)
    pression_Poisson(u,v,nx,ny,dx,dy,jmin,jmax,imax,imin)
    equilibre(u,v,u_star,v_star,p)
    save(record,t,u,v,p)
print("end Simulation !")


# plt.plot(record[3],record[0])
# plt.show()