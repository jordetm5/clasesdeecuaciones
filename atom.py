 
import numpy as np
import atomwrite


def PosicionInicial(N,L):
  Posicion= np.zeros((N,3), float)
  Nlat= int(N**(1./3.)+1.)
  r=L*(np.arange(Nlat,dtype=float)/Nlat-0.5) 
  i=0
  
  for x in r:
    for y in r:
      for z in r: 
          Posicion[i] =np.array([x,y,z],float)
          i+=1
        
          if i>=N:
            return Posicion  
  return Posicion
def inivel(n,t):
  vel=np.random.rand(n,3)
  return vel
def iniforces(pos,l):
  Forces=np.zeros_like(pos)
  return Forces
def Calforce(pos,vel,forces,l,dt):
  return 0

def update(v,p,d,dt,l,f):
  r=np.zeros([d,3],float)
  u=np.zeros([d,3],float)
  r[0]=p
  u[0]=v
  forces=iniforces(p,l)
  for j in range (0,3):
    for i in range (0, d-1):
      
      r[i+1][j] = r[i][j] + u[i][j] * dt + (Calforce(r[i][j],u[i][j],f,l,dt)* (dt**2)/2)
      u[i+1][j] = u[i][j] + ((Calforce(r[i+1][j],u[i+1][j],f,l,dt)+Calforce(r[i][j],u[i][j],f,l,dt))*dt/2)
      if r[i+1 ,j] > l or r[i+1 ,j] < 0 :
        r[i+1]=r[0]
     
  return r,u,forces

  
n=int(input("ingresa numero de moleculas: "))    
t=float(input("ingresa el tiempo: "))
dt=0.001
rho=0.005
L=(n/rho)**(1/3)
temp=1.0

d=int(t/dt)
p=PosicionInicial(n,L)
v=inivel(n,t)
f=iniforces(p,L)
l=np.arange(0.0,t,dt)
print(p)
print(v)

Pdb = atomwrite.pdbfile("animnew.pdb",L)
for i in range(0,d):
	s=np.zeros([n,3])
	for k in range(0,n):
  		x,y,z=update(v[k],p[k],d,dt,L,f)
		
		s[k]=x[i]
	
	Pdb.write(s)	
          


Pdb.close() 

