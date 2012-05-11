#
# source config.bash
# from homog import *
# (cell,mol) = domost()

from dolfin import *
import numpy as np
from params import *

# classes
from CellularDomain import *
from MolecularDomain import *
import scalar
import field 

class empty:pass


def u0_boundary(x, on_boundary):
  return on_boundary

## solv. homog cell
def solve_homogeneous_unit(problem,type="scalar"):

  ## debug mode
  debug=1
  if(debug==1):
    print "In debug mode - using stored values"
    d_eff = np.array([2.,2.,2.])
    problem.d_eff = d_eff
    V = FunctionSpace(problem.mesh,"CG",1)
    u = Function(V)
    u.vector()[:] = 1
    problem.x = u
    return 

  ## solve homog
  if(type=="scalar"):
    scalar.solveHomog(problem)
    D_eff = scalar.compute_eff_diff(problem)
    d_eff = np.array([D_eff,D_eff,D_eff])

  elif (type=="field"):
    field.solveHomog(problem)
    d_eff = field.compute_eff_diff(problem)

  else: 
    print "Not supported"
    quit()

  problem.d_eff = d_eff

## Coupled problem (need to check on neumann confs) 
# Here I'm just trying to get a general time-dep soln to work. NOthing is
# relevant to my problem yet 
def solve_homogenized_whole(cell,mol):

  # TODO using these for now, but not correct 
  # TODO probably want to pull values from 'x' used in homogeneous solution 
  u0 = 1.
  cell.bc = DirichletBC(cell.V, u0, u0_boundary)
  # TODO add in LCC instead
  dudn = -4  # assigning a flux on the entire boundary for now

  
  out  = File(problem.name+"_homogenized.pvd") 


  print "WARNING: overwriting anistropic D const. Needs to be fixed"
  print "WARNING: must adjust both unit cell and whole to enforce VecFunc basis"
  cell.d_eff = cell.d_eff[0]
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
  
  # for multiple bcs
  # for bc in bcs:
  #   bc.apply(A,b)

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

    print" TODO - note: ds needs to be limited to boundary between cell and microdomains [eq use ds(2)]"
    b[:] = 0.0
    # cell + molec coupling
    E = cell.d_eff * assemble( jcell * cell.v * ds, mesh=cell.mesh) 
    # outer cell boundary
    # TODO verify
    E = assemble(dcdn*cell.v*ds,mesh=cell.mesh)
    b[:] = E 
    

    b += M*x
    cell.bc.apply(A,b)
    solve(A,x,b,"gmres","default")
    #write

    cell.x.vector()[:] = x[:]


    # apparently each time we save, we get a new VTU for the time slice
    out << cell.x


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
  cell.mesh = UnitCube(8,8,8)    
  mol.mesh  = UnitCube(8,8,8)    
  
  solve_homogeneous_unit(cell)
  solve_homogeneous_unit(mol)
  
  compute_eff_diff(cell)
  compute_eff_diff(mol)

  return (cell,mol)

# test 2 
def Test():
  root = "/home/huskeypm/scratch/homog/mol/"

  # celular domain 
  prefix = "cell"
  meshFileOuter = root+prefix+"_mesh.xml.gz"
  subdomainFileOuter = root+prefix+"_subdomains.xml.gz"
  cellDomUnit = CellularDomain(meshFileOuter,subdomainFileOuter)
  cellDomUnit.Setup()
  cellDomUnit.AssignBC()
  solve_homogeneous_unit(cellDomUnit.problem,type="field")


  # molecular domain
  prefix = "mol"
  meshFileInner = root+prefix+"_mesh.xml.gz"
  subdomainFileInner = root+prefix+"_subdomains.xml.gz"
  molDomUnit = MolecularDomain(meshFileInner,subdomainFileInner)
  molDomUnit.Setup()
  molDomUnit.AssignBC()
  solve_homogeneous_unit(molDomUnit.problem,type="field")


  
  ## whole cell solutions 
  # solve on actual cell domain 
  # TODO Need to add
  print "Quit until whole cell mesh is loaded"
  prefix = "wholecell"
  meshFileOuter = root+prefix+"_mesh.xml.gz"
  subdomainFileOuter = root+prefix+"_subdomains.xml.gz"
  cellDomWhole = CellularDomain(meshFileOuter,subdomainFileOuter)
  cellDomWhole.Setup()
  cellDomWhole.AssignBC()
  cellDomWhole.problem.d_eff = cellDomUnit.problem.d_eff
  solve_homogenized_whole(cellDomWhile.problem,molDomUnit.problem)

  
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
  solve_homogenized_whole(cell,mol)

