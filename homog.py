#
# source config.bash
# from homog import *
# (cell,mol) = domost()

from dolfin import *
import numpy as np
from params import *

# classes
from CellularDomain import *
from CellularUnitDomain import *
from MolecularUnitDomain import *
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

def getb(problem):
  print "Replace w real flux"
  mesh = problem.mesh
  n = FacetNormal(mesh)
  flux = Constant(4.)
  #E = flux*dot(tetrahedron.n,v)*ds
  E = flux*dot(n,problem.v)*ds
  b  = assemble(E)
  return b
 
## Coupled problem (need to check on neumann confs)
# NOTE: assuming molecular domain is in steady state, so only solving
# cellular domain in time-dep fashion
def solve_homogenized_whole(wholecell,unitcell,unitmol,type="scalar"):
  # make Diff matrix
  Dii  = Constant((wholecell.d_eff[0],wholecell.d_eff[1],wholecell.d_eff[2]))
  Dij = diag(Dii)  # for now, but could be anisotropic


  # store output 
  out  = File(wholecell.name+"_homogenized.pvd")

  # Assembly of the K, M and A matrices
  wholecell.u = TrialFunction(wholecell.V)
  wholecell.v = TestFunction(wholecell.V)

  u,v = wholecell.u,wholecell.v
  a =inner(Dij * grad(u), grad(v))*dx
  K = assemble(a,mesh=wholecell.mesh)
  # TODO check on this
  M = assemble(inner(u,v)*dx,mesh=wholecell.mesh)

  u_n = Function(wholecell.V)
  A = K.copy()

  # TODO check on this 
  b = getb(wholecell)
  #print "Replace w real flux"
  #mesh = wholecell.mesh
  #n = FacetNormal(mesh)
  #flux = Constant(4.)
  ##E = flux*dot(tetrahedron.n,v)*ds
  #E = flux*dot(n,v)*ds
  #b  = assemble(E)

  # init cond
  x = u_n.vector()
  x[:] = 0.1 # initial conc  
  wholecell.x = Function(wholecell.V) # not sure why I did this 
  wholecell.x.vector()[:] = x[:]

  #b += M * x
  #quit()

  # for multiple bcs (though currently we do not have any dirichlet)
  for bc in wholecell.bcs:
     bc.apply(A,b)

  # time parms 
  dt =0.5
  t = 0.
  tMax = 1

  ## iterate
  while (t < tMax):

    ## TODO check that flux is correct
    # assume simple flux between compartments
    field.CalcConc(wholecell)
    print "WARNING: cmtd out"
    #field.CalcConc(unitmol)
    unitmol.conc=1

    # TODO - flux between unitmol and wholecell domain (where is the surface here, since this
    # concept seems tobe relevant only for unit cell??) 
    k = 1
    hack = 0.5
    jcell = k*(wholecell.conc - hack*unitmol.conc)
    print "wholecell: %f" % wholecell.conc
    print "unitmol: %f" % unitmol.conc
    jmol  = -jcell

    # TODO add scale factors (are these consistent)
    print "WARNING: cmtd out"
    #jcell*= unitcell.beta / unitcell.gamma
    #jmol *= unitmol.beta / unitmol.gamma
    
    ## TODO add corret time stepping, equations
    # solve cell
    # JOHAN
    print "t %f" % t
    t  += float(dt)
    A.assign(K)
    A *= parms.d*dt
    A += M

    #  TODOGoel: how/where is the surface defined here? It seems like the boundaries within
    # the unit cell are 'embedded' into the large cell description, so are no longer
    # boundaries
    print" TODO - note: ds needs to be limited to boundary between cell and microdomains [eq use ds(2)]"

    # cell + unitmolec coupling
    #print"FOR NOW IGNORING MOLEC CONTRIB"
    #E = cell.d_eff * assemble( jcell * cell.v * ds, mesh=cell.mesh)
    #b += E

    # outer cell boundary
    # TODO check on this 
    b1 = getb(wholecell)
    b=b1

    # add in mass matri
    b += M*x    

    #if(hasattr(wholecell,'bc')):
    #  wholecell.bc.apply(A,b)
    for bc in wholecell.bcs:
      bc.apply(A,b)

    ## solve and store 
    solve(A,x,b,"gmres","default")

    # store solution 
    wholecell.x.vector()[:] = x[:]

    # store results 
    out << wholecell.x


    # solv unitmol
    #F = unitmol.diff_eff * grad(unitmol.u) * grad(unitmol.v)
    #F += unitmol.diff_eff( junitmol * unitmol.v)
    #M = 0
    #solve(F==M,unitmol.u)
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
  cellDomUnit = CellularUnitDomain(meshFileOuter,subdomainFileOuter)
  cellDomUnit.Setup(type="field")
  cellDomUnit.AssignBC()
  solve_homogeneous_unit(cellDomUnit.problem,type="field")


  # molecular domain
  prefix = "mol"
  meshFileInner = root+prefix+"_mesh.xml.gz"
  subdomainFileInner = root+prefix+"_subdomains.xml.gz"
  molDomUnit = MolecularUnitDomain(meshFileInner,subdomainFileInner)
  molDomUnit.Setup(type="field")
  molDomUnit.AssignBC()
  solve_homogeneous_unit(molDomUnit.problem,type="field")



  ## whole cell solutions
  # solve on actual cell domain
  prefix = "wholecell"
  meshFileOuter = root+prefix+"_mesh.xml.gz"
  subdomainFileOuter = root+prefix+"_subdomains.xml.gz"
  cellDomWhole = CellularDomain(meshFileOuter,subdomainFileOuter)
  print "WARNING: using scalar for right now"
  cellDomWhole.Setup(type="field")
  # TODO: Need to replace, since using time-dep Neumann
  cellDomWhole.AssignBC()
  cellDomWhole.problem.d_eff = cellDomUnit.problem.d_eff
  solve_homogenized_whole(cellDomWhole.problem,cellDomUnit.problem,molDomUnit.problem)


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

  #GoelEx2p7()
  Test()
  quit()

  import sys
  if len(sys.argv) < 1:
      raise RuntimeError(msg)


  (cell,mol) = domost()
  solve_homogenized_whole(cell,mol)

