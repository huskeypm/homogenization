from dolfin import *
import numpy as np

boundsMin = np.zeros(3)
boundsMax = np.zeros(3)

EPS = 0.001 
def DefineBoundary(x,btype,on_boundary):
  if(not on_boundary):
    return 0 

  lr = ( np.abs(x[0]-boundsMin[0]) < EPS or np.abs(x[0]-boundsMax[0]) < EPS)
  tb = ( np.abs(x[1]-boundsMin[1]) < EPS or np.abs(x[1]-boundsMax[1]) < EPS)
  fb = ( np.abs(x[2]-boundsMin[2]) < EPS or np.abs(x[2]-boundsMax[2]) < EPS)
  obs = (not tb and not lr and not fb)

  if(btype=="lr"):
    return lr

  if(btype=="obs"):
    return obs

  if(btype=="tb"):
    return tb  

  if(btype=="fb"):
    return fb  

class XBoundary(SubDomain):
  def inside(self, x, on_boundary):
    #return (( x[0] < 0 or x[0] > 7) and on_boundary)
    return (DefineBoundary(x,"lr",on_boundary))

class YBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return (DefineBoundary(x,"tb", on_boundary))

class ZBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return (DefineBoundary(x,"fb", on_boundary))


class ObsBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return (DefineBoundary(x,"obs", on_boundary))



mesh = Mesh("3d.xml.gz") 

for i in np.arange(3):
  boundsMin[i] = np.min(mesh.coordinates()[:,i])
  boundsMax[i] = np.max(mesh.coordinates()[:,i])

print boundsMin
print boundsMax

# Test boundaries
test= 0
if (test==1):
  V = FunctionSpace(mesh,"CG",1)
  x = Function(V)
  bc_lr = DirichletBC(V,Constant(1),XBoundary())
  bc_lr.apply(x.vector())
  bc_obs = DirichletBC(V,Constant(20),ObsBoundary())
  bc_obs.apply(x.vector())
  bc_tb = DirichletBC(V,Constant(5),YBoundary())
  bc_tb.apply(x.vector())
  bc_tb = DirichletBC(V,Constant(7),ZBoundary())
  bc_tb.apply(x.vector())
  File("testbound.pvd") << x
  quit()
  


V = VectorFunctionSpace(mesh,"CG",1)

bcs = [] 
## e = [1,0,0]
bcs.append(DirichletBC(V.sub(0),Constant(0),XBoundary()))
# Neumann i believe is nalba w + delta = 0 on boundary, so natural
## e = [0,1,0]
bcs.append(DirichletBC(V.sub(1),Constant(0),YBoundary()))
# Neumann i believe is nalba w + delta = 0 on boundary, so natural
## e = [0,0,1]
bcs.append(DirichletBC(V.sub(2),Constant(0),ZBoundary()))
# Neumann i believe is nalba w + delta = 0 on boundary, so natural

Dii  = Constant((1.0,1.0,1.0))
Aij = diag(Dii) 
Delta = Identity( mesh.ufl_cell().geometric_dimension()) #

u = TrialFunction(V)
v = TestFunction(V)
gamer=1
if(gamer==1):
  form = inner(Aij*(grad(u) + Delta), grad(v))*dx(1)
else:
  form = inner(Aij*(grad(u) + Delta), grad(v))*dx
a = lhs(form)
L = rhs(form)

x = Function(V)
solve(a == L, x, bcs)
#plot(u,interactive=True)
File("3d.pvd") << x 

dim = mesh.ufl_cell().geometric_dimension()
omegas = np.zeros(dim)
for i in range(dim):
  grad_Xi_component = inner(grad(x[i]),Constant((1,0,0)))
  if(gamer==1):
    form = (grad_Xi_component + Constant(1)) * dx(1)
  else:
    form = (grad_Xi_component + Constant(1)) * dx
  integrand = assemble(form)
  omegas[i] = integrand

if(gamer==1):
  vol = assemble( Constant(1)*dx(1),mesh=mesh ) 
  surf = assemble( Constant(1)*ds,mesh=mesh ) 
else:
  vol = assemble( Constant(1)*dx,mesh=mesh ) 
  surf = assemble( Constant(1)*ds,mesh=mesh ) 
print omegas/vol


