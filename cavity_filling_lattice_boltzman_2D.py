import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib
import matplotlib.animation as animation
import matplotlib.patches as patches
from IPython import get_ipython
fs=20
plt.style.use('seaborn-dark')
plt.rc('xtick',labelsize=fs)
plt.rc('ytick',labelsize=fs)
plt.rc('text', usetex=True)
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.usetex'] = False


ca=np.array([[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]) #shéma D2Q9 donc 9 directions
w=[4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36] #racines des polynomes de Hermites (TP 0)
c0=1/np.sqrt(3) # facteur de normalisation des vitesses propres lorsque l'on conserve les moments jusqu'à l'ordre 2


# Maillage/Domaine

R=5 # rayon trou d'entrée

Lx=40*R # taille du domaine en fonction de R
Ly=20*R

nx=50*R # paramètre de maille, en fonction de R
ny=50*R

x0,y0=int(8*R),int(ny/2) # centre du cylindre

dx = R/nx # Taille de maille


# Vitesse du son réelle à 293K:
c0_real=np.sqrt(1.4*287*293) #racine de gamma*R*T pour gaz parfaits

dt=dx*c0/c0_real # Pas de temps
deltat=1 # delta de temps en unité réseau
M0=0.3 # Mach

U0=M0*c0 # vitesse du fluide initiale
print(U0)

print('nb de mailles : '+str(nx*ny))


def init(M0,Re):

    # Initialisation du domaine
    rho,uy=np.ones((nx,ny)),np.zeros((nx,ny))
    ux=U0*np.ones((nx,ny)) # le fluide est initalement à la vitesse U0 selon l'axe x
    nu=U0*2*R/Re #augmenter Re veut dire diminuer la viscosité
    taug=0.5+3*nu # on se sert de la formule du cours
    print('Valeur de taug : '+str(taug)+' \n')

    # on doit mettre à l'équilibre la fonction geq car le fluide n'est pas au repos
    geq=np.zeros((nx,ny,9))
    geq=eq(rho,ux,uy,geq,c0,ca)
    return geq,rho,ux,uy,taug,U0,nu

def eq(rho,ux,uy,geq,c0,ca):
    # Mise à jour de geq
    for i in range (0,9):
        # ca[i,0] donne le vecteur vitesse associé selon les polynomes de Hermite
        # développement polynomial de la fonction d'équilibre pour une quadrature de guass à l'ordre 2
        geq[:,:,i] = rho[:,:]*w[i]*(1.+ (ca[i][0]*ux + ca[i][1]*uy) / (c0**2) + ( (  ca[i][0]*ux + ca[i][1]*uy )**2) / (2.*c0**4) - (ux**2 + uy**2) / (2.*c0**2) )
    return geq

def propagate(g,gcoll):

    # Etape de propagation, identique au TP précédent

    g[:,:,0]=gcoll[:,:,0] #[0,0]

    g[1:,:,1]=gcoll[0:-1,:,1] #[1,0]
    g[:,1:,2]=gcoll[:,0:-1,2] #[0,1]
    g[0:-1,:,3]=gcoll[1:,:,3]#[-1,0]
    g[:,0:-1,4]=gcoll[:,1:,4]#[0,-1]

    g[1:,1:,5]=gcoll[0:-1,0:-1,5]#[1,1]
    g[0:-1,1:,6]=gcoll[1:,0:-1,6]#[-1,1]
    g[0:-1,:-1,7]=gcoll[1:,1:,7]#[-1,-1]
    g[1:,0:-1,8]=gcoll[0:-1,1:,8]#[1,-1]

    return(g)

def collide(gcoll,g,geq,taug):
    # Etape de collision
    gcoll[:,:,:]=g[:,:,:]-1/taug*(g[:,:,:]-geq[:,:,:])
    return(gcoll)

def macro(g,rho,ux,uy):

    # identique au TP précédent

    # calcul des variables macro
    rux,ruy,rho=np.zeros((nx,ny)),np.zeros((nx,ny)),np.zeros((nx,ny))
    for i in range (0,9):
        rho+=g[:,:,i]
        rux+=g[:,:,i]*ca[i,0]
        ruy+=g[:,:,i]*ca[i,1]
        #On est en quadrature d'ordre 3 pour la quadrature de Gauss, donc on a egalité des moments jusqu’à l’ordre 2 (modèle D2Q9)
    ux=rux/rho
    uy=ruy/rho

    return(rho,ux,uy)

def wall_noslip(gcoll,g,mask):
    # Paroi solide sur le cylindre: Bounce back avec frottement
    g[mask==1,0]=gcoll[mask==1,0]
    g[mask==1,1]=gcoll[mask==1,3]
    g[mask==1,2]=gcoll[mask==1,4]
    g[mask==1,3]=gcoll[mask==1,1]
    g[mask==1,4]=gcoll[mask==1,2]
    g[mask==1,5]=gcoll[mask==1,7]
    g[mask==1,6]=gcoll[mask==1,8]
    g[mask==1,7]=gcoll[mask==1,5]
    g[mask==1,8]=gcoll[mask==1,6]
    return(g)

def inflow(g,rho,ux,uy,M0,c0,ca,geq):

    # Conditions de vitesse au bord superieur
    uy[round(nx/2)-R:round(nx/2)+R,0]=-M0/np.sqrt(3)
    ux[round(nx/2)-R:round(nx/2)+R,0]=0


    # Macro: On impose ux et on en déduit rho_in:
    rho[0,:]=1/(1-ux[0,:])*(2*(g[0,:,6]+g[0,:,3]+g[0,:,7])+g[0,:,0]+g[0,:,2]+g[0,:,4])

    # Distributions: On impose les distributions inconnues avec le Bounce-Back hors equilibre

    #recalcul des équilibres pour rho_wall
    geq=eq(rho,ux,uy,geq,c0,ca) #recalcul des équilibres pour rho_inflow
    g[0,:,5]=g[0,:,7]+(geq[0,:,5]-geq[0,:,7])
    g[0,:,1]=g[0,:,3]+(geq[0,:,1]-geq[0,:,3])
    g[0,:,8]=g[0,:,6]+(geq[0,:,8]-geq[0,:,6])

    return(g,geq,ux,uy,rho)

def outflow(g,gcoll):
    # Conditions de gradient nul

    # Out
    g[-1,:,3]=gcoll[-2,:,3]
    g[-1,:,6]=gcoll[-2,:,6]
    g[-1,:,7]=gcoll[-2,:,7]
    # Top avec i=-1 et i-1=-2
    g[1:,-1,4]=gcoll[1:,-2,4]
    g[1:,-1,7]=gcoll[1:,-2,7]
    g[1:,-1,8]=gcoll[1:,-2,8]
    # Bottom pareil
    g[1:,0,2]=gcoll[1:,1,2]
    g[1:,0,5]=gcoll[1:,1,5]
    g[1:,0,6]=gcoll[1:,1,6]

    return(g)

# Marquage des conditions aux limites:



mask=np.zeros((nx,ny));

# on vectorise la création de l'obstacle
def obstacle(mask,R):
    mask[0,:]=1
    mask[-1,:]=1
    mask[:,-1]=1
    mask[0:round(nx/2)-R,0]=1
    mask[round(nx/2)+R:-1,0]=1

obstacle(mask,R)
plt.pcolormesh(mask,cmap='RdBu')
plt.colorbar()
plt.axis('equal')
plt.show()

# Nombre de Reynolds:
Re=300

# initialisation:
geq,rho,ux,uy,taug,U0,mu=init(M0,Re)
g,gcoll=geq.copy(),geq.copy()
nt=0
start=time.time()

nitération=100 # nombre d'itération souhaité
skip=nitération-50 # permet d'afficher à partir d'un certain nombre de tours

while (nt<nitération):

    nt+=1

    # boucle de LBM

    gcoll=collide(gcoll,g,geq,taug) # étape de collision

    g=propagate(g,gcoll) # étape de propagation
    g=wall_noslip(g,g,mask) # condition limite de bord d'obstacle
    g=outflow(gcoll,g) # condition limite de sortie Neumann
    g,geq,ux,uy,rho=inflow(g,rho,ux,uy,M0,c0,ca,geq) # condition limite entrée
    rho,ux,uy=macro(g,rho,ux,uy) # calcul des grandeurs macro

    geq=eq(rho,ux,uy,geq,c0,ca)


tcal=time.time()-start
if tcal>60:
    print(str(nt)+" itérations en "+str(int(tcal//60))+"m"+ str(int(tcal%60))+"s: Performances: "+str((1e-6)*nt*nx*ny/tcal)+"  MLUPS")
else :
    print(str(nt)+" itérations en "+str(tcal)+"s: Performances: "+str((1e-6)*nt*nx*ny/tcal)+"  MLUPS")


fig=plt.figure(figsize=(10,5)) # on veut afficher uy pour vérifier l'écoulement
plt.pcolormesh(uy/U0,cmap='RdBu')
plt.colorbar()
plt.xlabel(r'$Y/L$',fontsize=fs)
plt.ylabel(r'$X/L$',fontsize=fs)
plt.title(r'$Re='+str(Re)+'$',fontsize=fs)
plt.axis('equal')
plt.show()

fig=plt.figure(figsize=(10,5)) # on veut afficher le flot pour vérifier l'écoulement
x=np.linspace(0,nx,nx)
y=np.linspace(0,ny,ny)
plt.streamplot(x,y,(ux.T),uy.T)
plt.xlabel(r'$X/L$',fontsize=fs)
plt.ylabel(r'$Y/L$',fontsize=fs)
plt.title(r'$Re='+str(Re)+'$',fontsize=fs)
plt.axis('equal')
plt.show()



