#
# source config.bash
# from homog import *
# (cell,mol) = domost()

from dolfin import *
import numpy as np

# classes
from CellularDomain import *
from MolecularDomain import *

class empty:pass

parms = empty()
parms.d = 1.0 # diff. const

def u0_boundary(x, on_boundary):
  return on_boundary

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
def solve_homogeneous_unit(problem):
  from dolfin import nabla_grad as grad

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
  ## Plot solution and mesh
  #plot(u) plot(mesh)
  ## Dump solution to file in VTK format
  #file = File("poisson.pvd") file << u

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


## Coupled problem (need to check on neumann confs) 
# Here I'm just trying to get a general time-dep soln to work. NOthing is
# relevant to my problem yet 
def solve_homogenized(cell,mol):

  # TODO using these for now, but not correct 
  # TODO probably want to pull values from 'x' used in homogeneous solution 
  u0 = 1.
  cell.bc = DirichletBC(cell.V, u0, u0_boundary)

  


  # set up time dep form 
  
  # Assembly of the K, M and A matrices
  K = assemble(cell.d_eff * inner(grad(cell.u), grad(cell.v))*dx,mesh=cell.mesh)
  M = assemble(cell.u*cell.v*dx,mesh=cell.mesh)
  #E = assemble(-u*inner(a, grad(v))*dx,mesh=mesh)


# if(1):
#   print "Test, erase me"
#   V = FunctionSpace(cell.mesh, "CG", 1)
#   u = TrialFunction(V) 
#   v = TestFunction(V) 
#   K = assemble(cell.d_eff * inner(grad(u), grad(v))*dx,mesh=cell.mesh)
#   A = K.copy()
#   A.assign(K)
#   quit()
  

  u_n = Function(cell.V)
  A = K.copy()
  b = Vector(A.size(1))
  b[:]=0.0
  x = u_n.vector()
  x[:] = cell.x.vector()[:] # pass 'x' values from initial solution of homogeneous domains
  #cell.x = x
  #cell.u.vector()[:] = u0
  #mol.u.vector()[:] = 0.5*u0
  #mol.x = x
  

  dt =0.5 
  t = 0.
  tMax = 1
  while (t < tMax):
    print "t %f" % t

    t  += float(dt)

    ## TODO check that flux is correct
    # assume simple flux between compartments
    cell.conc = assemble( cell.x * dx,mesh=cell.mesh)     
    mol.conc = assemble( mol.x * dx,mesh=mol.mesh)     
    k = 1
    # TODO - need to understand how to get non-zero fluxes and why this matters
    hack = 0.5
    jcell = k*(cell.conc - hack*mol.conc)
    print "cell: %f" % cell.conc
    print "mol: %f" % mol.conc
    jmol  = -jcell

    # TODO add scale factors 


    ## TODO add corret time stepping, equations
    # solve cell
    # JOHAN 
    A.assign(K)
    A *= parms.d*dt
    A += M

    # TODO - note: ds needs to be limited to boundary between cell and microdomains [eq use ds(2)]
    b[:] = 0.0
    E = cell.d_eff * assemble( jcell * cell.v * ds, mesh=cell.mesh) 
    b[:] = E 
    

    b += M*x
    cell.bc.apply(A,b)
    solve(A,x,b,"gmres","default")
    #write

    cell.x.vector()[:] = x[:]


    # solv mol 
    #F = mol.diff_eff * grad(mol.u) * grad(mol.v) 
    #F += mol.diff_eff( jmol * mol.v) 
    #M = 0
    #solve(F==M,mol.u)
    #write


  quit()

def domost():
  parms.d = 1.
  cell = empty()
  mol  = empty()
  cell.name = "cell"
  mol.name = "mol"
  
  # TODO need to rescale size 
  cell.mesh = UnitCube(16,16,16)
  mol.mesh  = UnitCube(16,16,16) 
  
  solve_homogeneous_unit(cell)
  solve_homogeneous_unit(mol)
  
  compute_eff_diff(cell)
  compute_eff_diff(mol)

  return (cell,mol)

# test 2 
def Test():
  root = "/home/huskeypm/scratch/homog/mol/"

  ## unit cell solutions 
  # celular domain 
  prefix = "cell"
  meshFileOuter = root+prefix+"_mesh.xml.gz"
  subdomainFileOuter = root+prefix+"_subdomains.xml.gz"
  cellDomUnit = CellularDomain(meshFileOuter,subdomainFileOuter)
  cellDomUnit.Setup()
  cellDomUnit.AssignBC()
  solve_homogeneous_unit(cellDomUnit.problem)
  diff_eff = compute_eff_diff(cellDomUnit.problem)


  # molecular domain
  prefix = "mol"
  meshFileInner = root+prefix+"_mesh.xml.gz"
  subdomainFileInner = root+prefix+"_subdomains.xml.gz"
  molDomUnit = MolecularDomain(meshFileInner,subdomainFileInner)
  molDomUnit.Setup()
  molDomUnit.AssignBC()
  solve_homogeneous_unit(molDomUnit.problem)
  compute_eff_diff(molDomUnit.problem)

  
  ## whole cell solutions 
  # solve on actual cell domain 
  # TODO Need to add
  quit()
  prefix = "wholecell"
  meshFileOuter = root+prefix+"_mesh.xml.gz"
  subdomainFileOuter = root+prefix+"_subdomains.xml.gz"
  cellDomWhole = CellularDomain(meshFileOuter,subdomainFileOuter)
  cellDomWhole.Setup()
  cellDomWhole.AssignBC()
  cellDomWhole.problem.diff_eff = cellDom.problem.diff_eff
  solve_homogenized(cellDomWhile.problem,molDomUnit.problem)

  
  quit()

# example 2.7 in Goel paper 
def GoelEx2p7():
  ## micro 
  unitLength = np.array([0.1,0.1,1.0]) # length of microdomain in unit cell  
  diff = Constant(2.5) # isotropic  

  # our geometry needs to look something like the 'cross' in Fig 2 
  # creating this in Blender
  # mesh is crap, doesn't work 
  meshFile = "goel.xml.gz"
  problem = empty()
  problem.mesh = Mesh(meshFile)

  #mesh = UnitCube(6,6,6)
  #problem.name ="test"
  #problem.mesh = mesh 

  solve_homogeneous_unit(problem)
  quit()


  # note: authors take u0(x=[0.5,0.5,0,5]) = 0 to fix arb. const
  # i think I can treat this like an initial condition and set some 
  # location in u0 to be zero. 
  u0 = Constant(0)  # this doesn't do anything, but could use as templat ealter
  u_1 = interpolate(u0,V) # applies Expression u0 to FunctionSpace 
  # solve stuff
  # u_1.assign(u)
  

##
## MAIN
##

if __name__ == "__main__":
  msg="script.py <arg>"
  remap = "none"

  print "debug"
  #GoelEx2p7()
  Test()
  quit()

  import sys
  if len(sys.argv) < 1:
      raise RuntimeError(msg)
  

  (cell,mol) = domost()
  solve_homogenized(cell,mol)

