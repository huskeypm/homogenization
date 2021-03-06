"""
----------------------------------
Smolfin - a numerical solver of the Smoluchowski equation for interesting geometries
Copyright (C) 2012 Peter Kekenes-Huskey, Ph.D., huskeypm@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


----------------------------------------------------------------------------
"""
#
# Revisions
# 121009 fixed dimension error in DefineBoundary
#
from dolfin import *
import numpy as np

boundsMin = np.zeros(3)
boundsMax = np.zeros(3)
case="none"
meshFile="none"
dim=-1
gamer = 0

EPS = 0.001 
def DefineBoundary(x,btype,on_boundary):
  if(not on_boundary):
    return 0 

  # need this, since for some reason 'dim' is not global
  dim = np.shape(x)[0]


  lr = ( np.abs(x[0]-boundsMin[0]) < EPS or np.abs(x[0]-boundsMax[0]) < EPS)
  tb = ( np.abs(x[1]-boundsMin[1]) < EPS or np.abs(x[1]-boundsMax[1]) < EPS)
  if(dim==3):
    fb = ( np.abs(x[2]-boundsMin[2]) < EPS or np.abs(x[2]-boundsMax[2]) < EPS)
  else:
    fb=0
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

##
##
##



def doit(meshFile="none", case="none",gamer=0):
# my meshes
#case = "cheng"
#case = "auriault"
  if(case=="cheng"):
    # works 
    gamer=1
    mesh = Mesh("2d.xml.gz") 
    mesh = Mesh("3d.xml.gz") 
    #
    #??not sure what this was mesh = Mesh("/home/huskeypm/scratch/120816/test_mesh.xml.gz")
    for i in np.arange( mesh.ufl_cell().geometric_dimension() ):
      boundsMin[i] = np.min(mesh.coordinates()[:,i])
      boundsMax[i] = np.max(mesh.coordinates()[:,i])
  
  # aurialt test
  elif(case=="auriault"):
    gamer=0
    meshFile = "half_auriault.xml" # from Auriault paper
    meshFile = "half_auriault_3d.xml"
    mesh = Mesh(meshFile)
    boundsMin[:]=np.array([0,999,0])
    boundsMax[:]=np.array([5,999,5])
  
  # shorten geom
  elif(case=="shorten"):
    print "SHORTEN geometry"
    meshFile = "shorten.xml" # from Auriault paper
    #bounds[0,:] = np.array([0,6])
    #bounds[1,:] = np.array([0,3])    
    mesh = Mesh(meshFile)     
    gamer=0
    #
    for i in np.arange( mesh.ufl_cell().geometric_dimension() ):
      boundsMin[i] = np.min(mesh.coordinates()[:,i])
      boundsMax[i] = np.max(mesh.coordinates()[:,i])

    boundsMax[2]=1; boundsMin[2]=0
  
  else:
    mesh = Mesh(meshFile)     
    #gamer=0
    #
    for i in np.arange( mesh.ufl_cell().geometric_dimension() ):
      print "Marking direction %d" % i
      boundsMin[i] = np.min(mesh.coordinates()[:,i])
      boundsMax[i] = np.max(mesh.coordinates()[:,i])
  
  
  dim = mesh.ufl_cell().geometric_dimension()
  print "Using %d dim mesh" % dim
  print boundsMin
  print boundsMax
  
  # Test boundaries
  test= 0
  if (test==1):
    V = FunctionSpace(mesh,"CG",1)
    x = Function(V)
    #bc_lr = DirichletBC(V,Constant(1),XBoundary())
    #bc_lr.apply(x.vector())
    #bc_obs = DirichletBC(V,Constant(20),ObsBoundary())
    #bc_obs.apply(x.vector())
    #bc_tb = DirichletBC(V,Constant(5),YBoundary())
    #bc_tb.apply(x.vector())
    bc_fb = DirichletBC(V,Constant(7),ZBoundary())
    bc_fb.apply(x.vector())
    File("testbound.pvd") << x
    quit()
    
  
  
  V = VectorFunctionSpace(mesh,"CG",1)
  
  bcs = [] 
  if(case=="cheng"):
    ## e = [1,0,0]
    bcs.append(DirichletBC(V.sub(0),Constant(0),XBoundary()))
    # Neumann i believe is nalba w + delta = 0 on boundary, so natural
    ## e = [0,1,0]
    bcs.append(DirichletBC(V.sub(1),Constant(0),YBoundary()))
    # Neumann i believe is nalba w + delta = 0 on boundary, so natural
    ## e = [0,0,1]
    if(dim==3):
      bcs.append(DirichletBC(V.sub(2),Constant(0),ZBoundary()))
      # Neumann i believe is nalba w + delta = 0 on boundary, so natural
  
  elif(case=="auriault"):
    ## e = [1,0,0]
    bcs.append(DirichletBC(V.sub(0),Constant(0),XBoundary()))
    # Neumann i believe is nalba w + delta = 0 on boundary, so natural
   
    # none on Y direction  
   
    ## e = [0,0,1]
    if(dim==3):
      bcs.append(DirichletBC(V.sub(2),Constant(0),ZBoundary()))
    # Neumann i believe is nalba w + delta = 0 on boundary, so natural
  if(case=="shorten"):
    bcs.append(DirichletBC(V.sub(0),Constant(0),XBoundary()))
    bcs.append(DirichletBC(V.sub(1),Constant(0),YBoundary()))
    if(dim==3):
      bcs.append(DirichletBC(V.sub(2),Constant(0),ZBoundary()))
  else:
    bcs.append(DirichletBC(V.sub(0),Constant(0),XBoundary()))
    bcs.append(DirichletBC(V.sub(1),Constant(0),YBoundary()))
    if(dim==3):
      bcs.append(DirichletBC(V.sub(2),Constant(0),ZBoundary()))
    
  if(dim==2):
    Dii  = Constant((1.0,1.0))
  elif(dim==3):
    Dii  = Constant((1.0,1.0,1.0))
  Aij = diag(Dii) 
  #Delta = Identity( mesh.ufl_cell().geometric_dimension()) #
  Delta = Identity( dim )

  print "WARNING: CODE IS DEPREACATED"
  Vscalar = FunctionSpace(mesh,"CG",1)
  pmf = Function(Vscalar)
  pmf.vector()[:] = 1.
  Vmax = 999
  c = mesh.coordinates()
  idx=np.where(c[:,1] > 0.8)
  pmf.vector()[idx[0]]=Vmax
  idx=np.where(c[:,1] < 0.2)
  pmf.vector()[idx[0]]=Vmax
  
  u = TrialFunction(V)
  v = TestFunction(V)
  if(gamer==1):
    form = inner(Aij*pmf*(grad(u) + Delta), grad(v))*dx(1)
  else:
    form = inner(Aij*pmf*(grad(u) + Delta), grad(v))*dx
  a = lhs(form)
  L = rhs(form)
  
  x = Function(V)

  # NOT GOOD FOR MOST LARGE MESHES solve(a == L, x, bcs)
  # Use linear variational solver 
  # Equation is linear, since diffusion coefficient doesn't depend on u,
  lvproblem = LinearVariationalProblem(a,L, x, bcs=bcs)        
  solver = LinearVariationalSolver(lvproblem)
  solver.parameters["linear_solver"] = "gmres"
  solver.parameters["preconditioner"] = "ilu"
  solver.solve()

  #plot(u,interactive=True)
  name = "%dd.pvd" % dim
  File(name) << x 
  
  omegas = np.zeros(dim)
  for i in range(dim):
    v = [0,0,0]
    v[i] = 1
    if(dim==2):
      grad_Xi_component = pmf*inner(grad(x[i]),Constant((v[0],v[1]))) + Constant(1)
    elif(dim==3):
      grad_Xi_component = pmf*inner(grad(x[i]),Constant((v[0],v[1],v[2]))) + Constant(1)
    outname = "diff%d.pvd" % i
    Vscalar = FunctionSpace(mesh,"CG",1)
    File(outname)<<project(grad_Xi_component,V=Vscalar)
  
    if(gamer==1):
      form = grad_Xi_component * dx(1)
    else:
      form = grad_Xi_component * dx
    integrand = assemble(form)
    omegas[i] = integrand
  
  if(gamer==1):
    vol = assemble( Constant(1)*dx(1),mesh=mesh ) 
    surf = assemble( Constant(1)*ds,mesh=mesh ) 
  else:
    vol = assemble( Constant(1)*dx,mesh=mesh ) 
    surf = assemble( Constant(1)*ds,mesh=mesh ) 
  
  print "D*"
  diff=boundsMax-boundsMin
  if(dim==2):
    unitCellVol =   np.prod(diff[0:2])
  else:
    unitCellVol =   np.prod(diff[0:3])

  Ds = omegas/unitCellVol
  print Ds
  print "Dx"
  print Ds[0]

  print "Vol: %f SA: %f unitCellVol %f\n" % (vol, surf,unitCellVol)
  
  # for debug
  # Only relevant for cubes
  #print "Vol bounds  volfrac  2/3p 2p/2-p"
  #v = np.prod( (boundsMax-boundsMin)[0:dim])
  #vf = vol /v
  #print "%f %f %f %f" % (v,vf,2/3.*vf, 2*vf/(3-vf))

if __name__ == "__main__":
  msg="""
Purpose: 

Usage:
  homog3d.py <-case or -mesh>

  -case = default, auriault, cheng, shorten
Notes:
"""

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-case"):
      case=sys.argv[i+1]
    if(arg=="-mesh"):
      meshFile=sys.argv[i+1]
    if(arg=="-gamer"):
      gamer = 1



  #raise RuntimeError("error w rescaling diff. coef. by volume; fixed in homog.py but not here") 




  doit(case=case,meshFile=meshFile,gamer=gamer) 
