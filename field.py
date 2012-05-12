#
# Revisions
#       10.08.10 inception
#
#  source ~/sources/dolfin_smol/config.bash

from dolfin import *
from params import *


def boundary(x,on_boundary):
    return on_boundary

def CalcArea(problem):
  area = assemble(problem.x[0] * ds, mesh=problem.mesh)
  problem.area = area
  return area

# calculate concentration
def CalcConc(problem):
  problem.conc = assemble( problem.x[0] * dx,mesh=problem.mesh)
  problem.conc /= assemble( Constant(1)*dx,mesh=problem.mesh)
  

#
# Solve homogenized diff. eqn based on vector field 
#
## solv. homog cell
# See notetaker notes from 2012-04-19 (pg 11) 
# Using first-order representation from Higgins paper ()
# type - scalar is valid only for certain cases
def solveHomog(problem):
  # use class definition for BCs, etc
  if(hasattr(problem,'init')):
    print"Using class definition"
    print "IGNORING FOR NOW AND OVERWRITING"

  if(0):
    a =1 

  # load everything w defauit values 
  else:
    print "Overriding - WARNING: Should be in its own class"

    # Create mesh and define function space
    mesh = problem.mesh
    VV = VectorFunctionSpace(mesh, "Lagrange", 1)
    problem.VV = VV

    # Define boundary conditions
    # dirichlet 
    # TODO verify  
    u0 = Constant((1.,1.,1.))
    bc = DirichletBC(VV, u0, boundary)
    problem.bcs = [bc]

  # Define variational problem
  u = TrialFunction(VV)
  v = TestFunction(VV)
  
  ## LHS terms 
  # Diffusion constant
  #Dbulk = 2
  Dbulk = parms.d
  Dii  = Constant((Dbulk,Dbulk,Dbulk))
  Aij = diag(Dii)  # for now, but could be anisotropic
  
  # Identity matrix 
  Delta = Identity( mesh.ufl_cell().geometric_dimension()) #
  
  # LHS 
  form = inner(Aij*(grad(u) + Delta), grad(v))*dx
  
  # note: we are mixing linear and bilinear forms, so we need to split
  # these up for assembler 
  LHS = lhs(form)
  RHS = rhs(form)
  a = LHS
  
  ## RHS terms 
  print "WATNING: still confused about the BC here. Johan says ==0, Gary says we have complicated Jacobian term"
  n = FacetNormal(mesh)
  L = inner( n,v ) * ds 
  
  # add in RHS from earlier 
  L += RHS
  
  # Compute solution
  x = Function(VV)
  solve(a == L, x, bc)
  
  # Project solution to a continuous function space
  problem.x = x
  problem.up = project(x, V=VV)
  
  # save soln
  File(problem.name+"_unit.pvd") << problem.up

  return problem

def compute_eff_diff(problem):
  mesh = problem.mesh
  dim = mesh.ufl_cell().geometric_dimension()
  
  ## get omega
  import numpy as np
  omegas = np.zeros(dim)
  x = problem.x
  for i in range(dim):
    integrand = assemble( x[i]*dx)
    #print integrand 
    omegas[i] = integrand
  
  
  area = CalcArea()
  omegas /= 1  
  Deff = parms.d*omegas
  print Deff 
  
  return Deff

class empty:pass

def doit(fileIn):
  # Create mesh and define function space
  mesh = UnitCube(8,8,8)
  problem = empty()
  problem.mesh = mesh 
  solveHomog(problem)
 
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


