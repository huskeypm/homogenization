#
# source config.bash
# 

from dolfin import *
import numpy as np
from params import *
from dolfin import nabla_grad as grad


class empty:pass

# TODO - add boundary 
# might not be correct
def  area(mesh, boundary=-1):
  f = Constant(1)
  if(boundary==-1):
    A = f*ds
    a = assemble(A, mesh=mesh)
  else:
    A = f*ds(i)
    a = assemble(A, exterior_facet_domains=self.sub_domains, mesh=self.mesh)

  return a



## solv. homog cell
# See notetaker notes from 2012-04-19 (pg 11) 
# Using first-order representation from Higgins paper ()
# type - scalar is valid only for certain cases
def solveHomog(problem):
  # use class definition for BCs, etc
  if(hasattr(problem,'init')): 
    print"Using class definition"

  # load everything w defauit values 
  else:
    print "Overriding - WARNING: Should be in its own class"

    # Create mesh and define function space
    mesh = problem.mesh     
    V = FunctionSpace(mesh, "CG", 1)
    problem.V = V

    # Define boundary conditions
    # dirichlet 
    # TODO verify  
    u0 = Constant(1.)
    bc = DirichletBC(V, u0, u0_boundary)
    problem.bcs = [bc]

  # neumann cond
  # TODO verify 
  # JOHAN 
  dcdn = Constant(-1.0)

  # Define variational problem
  V = problem.V
  u = TrialFunction(V) 
  problem.u = u
  v = TestFunction(V) 
  problem.v = v
  a = parms.d * inner(grad(u), grad(v))*dx 

  f = Constant(0.0) 
  L = f*v*dx - dcdn*v*ds
  # Compute solution
  x = Function(V) 
  solve(a == L, x, problem.bcs)
  problem.x = x

  # save soln
  File(problem.name+"_unit.pvd") << problem.x
  

## compute effective diff const 
def compute_eff_diff(problem):
  from dolfin import assemble
  # TODO need to computre area 
  #gamma = assemble(ds)    

  gamma = area(problem.mesh,boundary=-1)

  # TODO need to figure out how to represent ident
  #I = Identity(problem.V.cell().d)

  # TODO FIXME     
  #I = delta
  # JOHAN 
  #integral = assemble( (grad(problem.u) + I ) * dx ) 
  print "not correct"
  integral = 1
  omegac = (1/gamma) * integral

  diffeff = omegac * parms.d
  problem.d_eff = diffeff 

  return diffeff



