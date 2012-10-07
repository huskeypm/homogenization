#
# Revisions
# 121004 expanded functionality to handle 'auriault' geometry from 97 paper
#        results seem reasonable 
#
from dolfin import *
import numpy as np



bounds = np.zeros([2,2])

EPS = 0.001 
def DefineBoundary(x,btype,on_boundary):
  if(not on_boundary):
    return 0 

  #lr = ( x[0] < 0 or x[0] > 7) 
  #print bounds  
  lr = ( x[0] < bounds[0,0]+EPS or x[0] > bounds[0,1]-EPS) 
  tb = (not lr and ( x[1] < bounds[1,0]+EPS or x[1] > bounds[1,1]-EPS) )
  obs = (not tb and not lr)

  if(btype=="lr"):
    return lr

  if(btype=="obs"):
    return obs

  if(btype=="tb"):
    return tb  

class XBoundary(SubDomain):
  def inside(self, x, on_boundary):
    #return (( x[0] < 0 or x[0] > 7) and on_boundary)
    return (DefineBoundary(x,"lr",on_boundary))

class YBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return (DefineBoundary(x,"tb", on_boundary))

class ObsBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return (DefineBoundary(x,"obs", on_boundary))


cheng=0
auriault=1
case ="shorten" # cheng, auriault
if(case=="cheng"):
  print "CHENG geometry"
  meshFile = "2d.xml" # original from Cheng paper, I think 
  bounds[0,:] = np.array([0,7])
  bounds[0,:] = np.array([1,5])
elif(case=="auriault"):
  print "AURIAULT geometry"
  meshFile = "half_auriault.xml" # from Auriault paper
  bounds[0,:] = np.array([0,5])
  bounds[1,:] = np.array([999,999])
elif(case=="shorten"):
  print "SHORTEN geometry"
  meshFile = "shorten.xml" # from Auriault paper
  #bounds[0,:] = np.array([0,6])
  #bounds[1,:] = np.array([0,3])    
  mesh = Mesh(meshFile)     
  #
  for i in np.arange(3):
    boundsMin[i] = np.min(mesh.coordinates()[:,i])
    boundsMax[i] = np.max(mesh.coordinates()[:,i])



mesh = Mesh(meshFile) 

# Test boundaries
test= 0
if (test==1):
  V = FunctionSpace(mesh,"CG",1)
  x = Function(V)
  bc_lr = DirichletBC(V,Constant(1),XBoundary())
  bc_lr.apply(x.vector())
  bc_obs = DirichletBC(V,Constant(10),ObsBoundary())
  bc_obs.apply(x.vector())
  bc_tb = DirichletBC(V,Constant(5),YBoundary())
  bc_tb.apply(x.vector())
  File("testbound.pvd") << x
  quit()
  


V = VectorFunctionSpace(mesh,"CG",1)

bcs = [] 
if(case=="cheng"):
  ## e = [1,0]
  bcs.append(DirichletBC(V.sub(0),Constant(0),XBoundary()))
  # Neumann i believe is nalba w + delta = 0 on boundary, so natural
  ## e = [0,1]
  bcs.append(DirichletBC(V.sub(1),Constant(0),YBoundary()))
  # Neumann i believe is nalba w + delta = 0 on boundary, so natural
elif(case=="auriault"):
  # left/right periodic 
  bcs.append(DirichletBC(V.sub(0),Constant(0),XBoundary()))
elif(case=="shorten"):
  # left/right periodic 
  bcs.append(DirichletBC(V.sub(0),Constant(0),XBoundary()))
  # Neumann i believe is nalba w + delta = 0 on boundary, so natural
  # top/bottom periodic 
  bcs.append(DirichletBC(V.sub(1),Constant(0),YBoundary()))

Dii  = Constant((1.0,1.0))
Aij = diag(Dii) 
Delta = Identity( mesh.ufl_cell().geometric_dimension()) #

u = TrialFunction(V)
v = TestFunction(V)
form = inner(Aij*(grad(u) + Delta), grad(v))*dx
a = lhs(form)
L = rhs(form)

x = Function(V)
solve(a == L, x, bcs)
#plot(u,interactive=True)
File("2d.pvd") << x 

dim = 2 
omegas = np.zeros(dim)
for i in range(dim):
  v = [0,0]
  v[i] = 1
  grad_Xi_component = inner(grad(x[i]),Constant((v[0],v[1]))) + Constant(1)
  outname = "diff2d%d.pvd" % i
  Vscalar = FunctionSpace(mesh,"CG",1)
  File(outname)<<project(grad_Xi_component,V=Vscalar)

  form = grad_Xi_component * dx
  integrand = assemble(form)
  omegas[i] = integrand

vol = assemble( Constant(1)*dx,mesh=mesh ) 
print "D*"
print omegas/vol


