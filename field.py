#
# Revisions
#       10.08.10 inception
#
#  source ~/sources/dolfin_smol/config.bash

from dolfin import *
from params import *
from util import *



EPS = 1.e-10
################/////////////////


# calculate concentration
def CalcConc(domain):
  problem = domain.problem
  if(domain.gamer==0):
    problem.conc = assemble( problem.x[0] * dx,mesh=problem.mesh)
  else:
    problem.conc = assemble( problem.x[0] * dx(1),mesh=problem.mesh)

  #problem.conc /= assemble( Constant(1)*dx,mesh=problem.mesh)
  problem.conc /= problem.volume
  

#
# Solve homogenized diff. eqn based on vector field 
#
## solv. homog cell
# See notetaker notes from 2012-04-19 (pg 11) 
# Using first-order representation from Higgins paper ()
# type - scalar is valid only for certain cases
def solveHomog(domain):
  # mesh
  problem = domain.problem
  mesh = problem.mesh
  V = problem.V

  ################## /////////////////////
  prob = problem 
  problem.meshType = "gamer" # gamer
  print "NOT CORRECT REMOV EME"
  utilObj = util(problem)
  utilObj.GeometryInitializations(mesh)
  utilObj.DefinePBCMappings(mesh)

  Vv = VectorFunctionSpace(mesh, "CG", 1)
  class LeftRightBoundary(SubDomain):
      def inside(self, x, on_boundary):
          # find v1x
          return tuple(x) in problem.targetsx
      
      
      # field component 1
      def map(self, x, y):
          y[:] = problem.vert_mapx.get(tuple(x), x)
      
  class BackFrontDomain(SubDomain):
      def inside(self, x, on_boundary):
          # find v1x
          return tuple(x) in problem.targetsy
    
    
      # field component 1
      def map(self, x, y): 
          y[:] = problem.vert_mapy.get(tuple(x), x)
    
      
  class TopBottomDomain(SubDomain):
      def inside(self, x, on_boundary):
          # find v1x
          return tuple(x) in problem.targetsz
    
    
      # field component 1
      def map(self, x, y):
          y[:] = problem.vert_mapz.get(tuple(x), x)
    
  # Sub domain for Dirichlet boundary condition
  class CenterDomain(SubDomain):
    
      def inside(self, x, in_boundary):
          return all(near(x[i], problem.center_coord[i], EPS) for i in range(problem.nDims))


  # Create Dirichlet boundary condition
  bcs = []
  fixed_center = DirichletBC(Vv, Constant((0,0,0)), CenterDomain(), "pointwise")
  bcs.append(fixed_center)

  bc1 = PeriodicBC(Vv.sub(0), LeftRightBoundary())
  bcs.append(bc1)
  bc2 = PeriodicBC(Vv.sub(1), BackFrontDomain())
  bcs.append(bc2)
  bc3 = PeriodicBC(Vv.sub(2), TopBottomDomain())
  bcs.append(bc3)

  testBC = 0
  if(testBC):
    print "TEsting BCs"
    bc1 =DirichletBC(Vv, Constant((1,1,1)),LeftRightBoundary())
    bc2 =DirichletBC(Vv, Constant((1,1,1)),BackFrontDomain())
    bc3 =DirichletBC(Vv, Constant((1,1,1)),TopBottomDomain())
    bcs = []
    bcs.append(fixed_center)
    #bcs.append(bc1)
    #bcs.append(bc2)
    bcs.append(bc3)

  # Define variational problem
  u = TrialFunction(Vv)
  v = TestFunction(Vv)
  Dbulk = 2.
  Dii  = Constant((Dbulk,Dbulk,Dbulk))
  Aij = diag(Dii)
  Delta = Identity( mesh.ufl_cell().geometric_dimension()) #
  if(problem.meshType == "fenics"):
    form = inner(Aij*(grad(u) + Delta), grad(v))*dx
  elif(problem.meshType == "gamer"):
    form = inner(Aij*(grad(u) + Delta), grad(v))*dx(1)
  a = lhs(form)
  L = rhs(form)

  # Compute solution
  u = Function(Vv)
  solve(a == L, u, bcs)

  #plot(u)
  #plot(mesh)
  #interactive()

  # Save solution to file
  file = File("periodic.pvd")
  file << u




  print "NOT CORRECT REMOV EME"
  quit()
  ###################//////////////////////////

  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  
  ## LHS terms 
  # Diffusion constant
  Dbulk = parms.d
  print "OVERRIDING Dbulk"
  Dbulk = 2.
  Dii  = Constant((Dbulk,Dbulk,Dbulk))
  Aij = diag(Dii)  # for now, but could be anisotropic
  
  # Identity matrix 
  Delta = Identity( mesh.ufl_cell().geometric_dimension()) #
  
  # LHS 
  if(domain.gamer==0):
  #form = inner(Aij*(grad(u) + Delta), grad(v))*dx
    form = inner(Aij*(grad(u) + Delta), grad(v))*dx
  else:
    form = inner(Aij*(grad(u) + Delta), grad(v))*dx(1) 
  
  # note: we are mixing linear and bilinear forms, so we need to split
  # these up for assembler 
  LHS = lhs(form)
  RHS = rhs(form)
  a = LHS
  
  ## RHS terms 
  #print "WATNING: still confused about the BC here. Johan says ==0, Gary says we have complicated Jacobian term"
  #n = FacetNormal(mesh)
  #L = inner( n,v ) * ds 
  # Based on Johan's May 14 email, because (d X + I) .n = 0, we do not apply the BC. 
  
  # add in RHS from earlier 
  L = RHS

  
  # Compute solution
  x = Function(V)
  solve(a == L, x, problem.bcs)

  # not needewd
  file = File("periodic.pvd")
  file << x

  
  # Project solution to a continuous function space
  problem.x = x
  problem.up = project(x, V=V)
  
  # save soln
  File(problem.name+"_unit.pvd") << problem.up

  return problem

def compute_eff_diff(domain):
  problem = domain.problem
  mesh = problem.mesh
  dim = mesh.ufl_cell().geometric_dimension()
  
  ## get omega
  # treating each component independtly, following Goel's example in sec 2.7 
  import numpy as np
  omegas = np.zeros(dim)
  x = problem.x
  for i in range(dim):
    #form = (inner(grad(x[i]),Constant((1,0,0)))+Constant(1)) * dx # works 
    grad_Xi_component = inner(grad(x[i]),Constant((1,0,0)))

    if (domain.gamer==0):
      form = (grad_Xi_component + Constant(1)) * dx 
    else:
      form = (grad_Xi_component + Constant(1)) * dx(1)

    integrand = assemble(form)
    omegas[i] = integrand
  
  
  omegas /= problem.gamma
  d_eff = parms.d*omegas
  print "d_eff:"
  print d_eff

  print "rewighting by unit cell vol"
  d_eff /= problem.volUnitCell
  print d_eff

  # store 
  problem.d_eff = d_eff
  
  return d_eff

class empty:pass

def doit(fileIn):
  # Create mesh and define function space
  defaultDom = DefaultUnitDomain()
  mesh = UnitCube(8,8,8)
  problem = empty()
  problem.mesh = mesh 
  solveHomog(defaultDom,type="field")
 
import sys

if __name__ == "__main__":
  msg="Purpose: Runs diff. eqn with field as test function "
  msg=msg+"Usage: "
  msg=msg+"script.py <arg>"
  msg=msg+"Notes:"
  remap = "none"


  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"




  doit(fileIn)


