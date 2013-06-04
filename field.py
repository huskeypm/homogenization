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
#       10.08.10 inception
#
#  source ~/sources/dolfin_smol/config.bash

from dolfin import *
from params import *
from util import *
import numpy as np
import numpy as np



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
def solveHomog(domain,smolMode="false"):
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
  nDims = problem.nDims
  Dii  = Constant(Dbulk*np.ones(nDims))
  Aij = diag(Dii)  # for now, but could be anisotropic
  Atilde = Aij 

  # electro 
  if(smolMode!="false"):
    if(domain.gamer==1):
      raise RuntimeError("project() does not work with Gamer meshes. Try removing 'domain' info from mesh and rerunning without gamer tag")
    print "Adding in electrostatic component" 
    Vscalar = FunctionSpace(mesh,"CG",1)
    intfact = Function(Vscalar)
    intfact.vector()[:]    =    np.exp(-parms.beta * problem.pmf.vector()[:])
    Atilde = intfact * Aij
    File("distro.pvd") << intfact
  
  # Identity matrix 
  Delta = Identity( mesh.ufl_cell().geometric_dimension()) #
 
  # LHS 
  #print "Gamer", domain.gamer
  if(domain.gamer==0):
    form = inner(Atilde*(grad(u) + Delta), grad(v))*dx
  else:
    form = inner(Atilde*(grad(u) + Delta), grad(v))*dx(1) 
  
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
  #solve(a == L, x, problem.bcs)
  lvproblem = LinearVariationalProblem(a,L, x, bcs=problem.bcs)
  solver = LinearVariationalSolver(lvproblem)
  solver.parameters["linear_solver"] = "gmres"
  solver.parameters["preconditioner"] = "ilu"
  print "Using amg preconditioner instead of ilu"
  solver.parameters["preconditioner"] = "amg"
  solver.solve()

  problem.x = x
  problem.up = Function(problem.V)   


  if(smolMode!="false"):
    debugin=1
    if(debugin):

      File("xi.pvd") << project(problem.x[0],V=Vscalar)
      File("int.pvd") << project(intfact,V=Vscalar)
      File("prod.pvd") << project(problem.x[0]*intfact,V=Vscalar)
      File("pmf.pvd") << project(problem.pmf,V=Vscalar)
      p = problem.x[0]*intfact
      grad_Xi_component = p.dx(1)
      File("deriv.pvd") << project(grad_Xi_component, V=Vscalar)
      print "testing"

    print "WARNING: need to verify that electrostatic part applied correctly for homogeniztion"
    for i in range(problem.nDims):
      id = "%d" % i
      n = problem.x[i] * intfact
      n = problem.x[i]
      temp = project(n,V=Vscalar)
      fileName = problem.name+"_pmfprojected"+id+".pvd"
      File(fileName) << temp 
   
      # DEBUG 
      #test = project(x[i],V=Vscalar)
      #fileName = problem.name+"_comp"+id+".pvd"
      #File(fileName) << test
      ## NO WORK fileName = problem.name+"_test2.pvd"
      ##File(fileName) << x[0]

      # apply 
      #problem.x[i].vector()[:] = n.vector()[:]
      #ci = project(problem.x[i]*problem.pmf,V=Vscalar)
      #l = len(ci.vector())
      #print i,l
      #print len(problem.up.vector())
      #problem.up.vector()[i*l:(i+1)*l] = ci.vector()

    print "WARNING: may have to redo 'loop' over x[i] to get projection right"
    problem.up = project(problem.x * intfact, V=V)

    #problem.up[1] = project(problem.x[1],V=Vscalar)


    #print "ERROR: I CHANGED THIS TO GIVE WRONG ANSWER"
    #x.vector()[:] = x.vector()[:] * 1.
    #problem.up = x

  # leave
  else:
    problem.up.vector()[:] = problem.x.vector()[:]
    ## Project solution to a continuous function space
    # WARNING: this projection doesn't seem to give me meaningful output (e.q solns == 0)
    #problem.up = project(x, V=V)
  

  
  # save soln
  #File(problem.name+"_unit.pvd") << problem.up

  # save unprojected soln instead 
  fileName = problem.outpath + problem.name+"_unit.pvd"
  print "Writing ",fileName
  #File(fileName) <<  problem.x   
  File(fileName) <<  problem.up


  return problem

def compute_eff_diff(domain):
  problem = domain.problem
  mesh = problem.mesh
  dim = mesh.ufl_cell().geometric_dimension()


  if(domain.gamer==1):
    dx_int = dx(1)
  else:
    dx_int = dx

  Vscalar = FunctionSpace(mesh,"CG",1)
  us = TrialFunction(Vscalar)
  vs = TestFunction(Vscalar)

  ## get omega
  # treating each component independtly, following Goel's example in sec 2.7 
  import numpy as np
  omegas = np.zeros(dim)
  x = problem.up
  # Believe it is correct now print "WARNING: verify accuracy of this calculation"
  # I iterate only over the diagonal, since my original diffusion constant is diagonal 
  for i in range(dim):
    #form = (inner(grad(x[i]),Constant((1,0,0)))+Constant(1)) * dx # works 
    # Here  we suppose to be extracting the ith derivative of x[i]
    #v = [0,0,0]
    #v[i] = 1
    #grad_Xi_component = inner(grad(x[i]),Constant((v[0],v[1],v[2]))) + Constant(1)

    # test 
    #outname = "diff%d.pvd" % i
    #File(outname)<<project(grad_Xi_component,V=Vscalar)

    #if (domain.gamer==0):
    #  form = grad_Xi_component * dx 
    #else:
    #  form = grad_Xi_component * dx(1)
#
#    integrand = assemble(form)
#    omegas[i] = integrand

    # JOHAN 
    grad_Xi_component = x[i].dx(i)+Constant(1)
    outname = "diff%d.pvd" % i

    D_eff_project = Function(Vscalar)
   
    solve(us*vs*dx_int==grad_Xi_component*vs*dx_int, D_eff_project)
    #File(outname)<<D_eff_project
  
    form = grad_Xi_component * dx_int
  
    omegas[i] = assemble(form)
  
  vol = assemble( Constant(1)*dx_int, mesh=mesh )
  

  print "omegasO"
  print omegas
  #print omegas/vol
  #problem.omv = omegas/vol
  #problem.domv = parms.d*omegas/vol
  #problem.d_eff = problem.domv 
  #print "WARNING: circumventing use of volUnitCell etc for determining deff"
  #return problem.d_eff 
  #print "vol", vol
  
  
  print "WARNING: what is gamma?" 
  #omegas /= problem.gamma
  print "omegas"
  print omegas
  d_eff = parms.d*omegas

  print "Reweighting by unit cell vol"
  d_eff /= problem.volUnitCell
  if(dim==3):
    print "d_eff= [%4.2f,%4.2f,%4.2f] for d=%4.2f"%\
      (d_eff[0],d_eff[1],d_eff[2],parms.d)
  else: 
    print "d_eff= [%4.2f,%4.2f] for d=%4.2f"%\
      (d_eff[0],d_eff[1],parms.d)
  print "problem.volUnitCell", problem.volUnitCell

  # normalize
  #nd_eff= d_eff/ np.linalg.norm(d_eff)
  #print "Deff (normalized) ", nd_eff

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


