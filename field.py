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
# smol - add in smol. approach 
def solveHomog(domain,smol="false"):
  # mesh
  problem = domain.problem
  mesh = problem.mesh
  V = problem.V

  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  
  ## LHS terms 
  # Diffusion constant
  Dbulk = parms.d
  Dii  = Constant((Dbulk,Dbulk,Dbulk))
  Aij = diag(Dii)  # for now, but could be anisotropic

  if(smol!="false"):
    print "Need to add smol"
    quit()
  
  # Identity matrix 
  Delta = Identity( mesh.ufl_cell().geometric_dimension()) #
 
  # LHS 
  #print "Gamer", domain.gamer
  if(domain.gamer==0):
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
  problem.x = x

  
  ## Project solution to a continuous function space
  # WARNING: this projection doesn't seem to give me meaningful output (e.q solns == 0)
  #problem.up = project(x, V=V)
  
  # save soln
  #File(problem.name+"_unit.pvd") << problem.up

  # save unprojected soln instead 
  fileName = problem.name+"_unit.pvd"
  print "Writing ",fileName
  File(fileName) <<  x


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

  print "Reweighting by unit cell vol"
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


