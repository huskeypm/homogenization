#
# source config.bash
# from homog import *
# (cell,mol) = domost()

from dolfin import *
import numpy as np
from params import *

# classes
from CellularDomain import *
from DefaultUnitDomain import *
from CellularUnitDomain import *
from MolecularUnitDomain import *
import scalar
import field



class empty:pass


## solv. homog cell
def solve_homogeneous_unit(problem,type="scalar",debug=0):

  ## debug mode
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
  # using scalar fields 
  if(type=="scalar"):
    print "I cannot guarantee this is correct..."
    scalar.solveHomog(problem)
    D_eff = scalar.compute_eff_diff(problem)
    d_eff = np.array([D_eff,D_eff,D_eff])

  # using vector fields 
  elif (type=="field"):
    field.solveHomog(problem)
    d_eff = field.compute_eff_diff(problem)

  else:
    print "Not supported"
    quit()

  problem.d_eff = d_eff

# add in (time-dependent) flux condition
def getb(problem,t=0):
  print "Replace w real flux"
  mesh = problem.mesh
  n = FacetNormal(mesh)


  ## NOTE: I think it makes sense to keep this as simple as possible (e.g. no complicated fluxes or markers) 
  ## since mostly likely what I want is already in subcell. 
  flux = problem.dudn # hopefull the expression is getting updated with each time step 
  E = flux*dot(n,problem.v)*ds
  b  = assemble(E)
  return b
 
## Coupled problem (need to check on neumann confs)
# NOTE: assuming molecular domain is in steady state, so only solving
# cellular domain in time-dep fashion
def solve_homogenized_whole(wholecell,unitcell,unitmol,type="scalar",debug=0):


  # store output 

  ## molecular domain (steady state)
  # make Diff matrix
  print "NEED TO SOLVE MOL"
##  print "Solving molecular steady state..."
##  Dii  = Constant((wholecell.d_eff[0],wholecell.d_eff[1],wholecell.d_eff[2]))
##  Dij = diag(Dii)  # for now, but could be anisotropic
##  u,v = TrialFunction(unitmol.V), TestFunction(unitmol.V)
##  A = inner(Dij*grad(u), grad(v))*dx
##  L = 0
##  x = Function(unitmol.V)
##  
##  bcs = unitmol.bcs # use same Dirichlet cond. from unitcell prob.
##  # no flux bc on surface, so we don't explicitly add it here
##  solve(A==L,x,bcs) 
##  File(unitmol.name+"_homogenized.pvd") << x


  ## Wholecell
  print "Solving time-dependent cellular problem"
  # make Diff matrix
  Dii  = Constant((wholecell.d_eff[0],wholecell.d_eff[1],wholecell.d_eff[2]))
  Dij = diag(Dii)  # for now, but could be anisotropic
  # Assembly of the K, M and A matrices
  u,v = TrialFunction(wholecell.V), TestFunction(wholecell.V)
  wholecell.u = u
  wholecell.v = v
  a =inner(Dij * grad(u), grad(v))*dx
  K = assemble(a,mesh=wholecell.mesh)
  A = K.copy()
  # TODO check on this
  M = assemble(inner(u,v)*dx,mesh=wholecell.mesh)

  # TODO check on this 
  b = getb(wholecell)
  # for multiple bcs (though currently we do not have any dirichlet)
  for bc in wholecell.bcs:
     bc.apply(A,b)


  # init cond
  u_n = Function(wholecell.V)
  x = u_n.vector()
  x[:] = 0.1 # initial conc  
  wholecell.x = Function(wholecell.V) # not sure why I did this 
  wholecell.x.vector()[:] = x[:]

  out  = File(wholecell.name+"_homogenized.pvd")



  # time parms 
  dt =0.005
  t = 0.
  if(debug ==0):
    tStep = 10
  else:
    tStep= 2
  
  tMax = tStep * dt
  

  ## iterate
  while (t < tMax):

    ## TODO check that flux is correct
    # assume simple flux between compartments
    field.CalcConc(wholecell)
    if(debug==0):
      field.CalcConc(unitmol)
    else:
      unitmol.conc=0.2

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

##def domost():
##  parms.d = 1.
##  cell = empty()
##  mol  = empty()
##  cell.name = "cell"
##  mol.name = "mol"
##
##  # TODO need to rescale size
##  cell.mesh = UnitCube(8,8,8)
##  mol.mesh  = UnitCube(8,8,8)
##
##  solve_homogeneous_unit(cell)
##  solve_homogeneous_unit(mol)
##
##  compute_eff_diff(cell)
##  compute_eff_diff(mol)
##
##  return (cell,mol)

def Debug():
  cellDomUnit = DefaultUnitDomain()
  cellDomUnit.Setup(type="field")
  cellDomUnit.AssignBC()

  solve_homogeneous_unit(cellDomUnit.problem,type="field")

def Debug2():
  1 

# test 2
def SolveHom(root="./",debug=0):

  ## Setup
  root = "/home/huskeypm/scratch/homog/mol/"

  # celular domain
  prefix = "cell"
  meshFileOuter = root+prefix+"_mesh.xml.gz"
  subdomainFileOuter = root+prefix+"_subdomains.xml.gz"
  cellDomUnit = CellularUnitDomain(meshFileOuter,subdomainFileOuter)
  cellDomUnit.Setup(type="field")
  cellDomUnit.AssignBC()

  # molecular domain
  prefix = "mol"
  meshFileInner = root+prefix+"_mesh.xml.gz"
  subdomainFileInner = root+prefix+"_subdomains.xml.gz"
  molDomUnit = MolecularUnitDomain(meshFileInner,subdomainFileInner)
  molDomUnit.Setup(type="field")
  molDomUnit.AssignBC()

  # get fractional volume 
  cellProblem = cellDomUnit.problem
  molProblem = molDomUnit.problem
  volUnitCell = cellProblem.volume + molProblem.volume 
  cellProblem.gamma = cellProblem.volume/volUnitCell
  cellProblem.volUnitCell = volUnitCell
  print "cell frac vol: %f" % cellProblem.gamma
  molProblem.gamma = molProblem.volume/volUnitCell
  molProblem.volUnitCell = volUnitCell
  print "mol frac vol: %f" % molProblem.gamma


  ## Solve unit cell problems  
  print "Solving cellular unit cell using %s" % meshFileOuter
  solve_homogeneous_unit(cellDomUnit.problem,type="field",debug=debug)
  print "Solving molecular unit cell using %s"% meshFileInner
  solve_homogeneous_unit(molDomUnit.problem,type="field",debug=debug)


  ## whole cell solutions
  # solve on actual cell domain
  prefix = "multi_clustered"
  meshFileCellular= root+prefix+"_mesh.xml.gz"
  subdomainFileCellular= root+prefix+"_subdomains.xml.gz"
  if(debug==0):
    cellDomWhole = CellularDomain(meshFileCellular,subdomainFileCellular)
  else:
    cellDomWhole = DefaultUnitDomain()

  cellDomWhole.Setup(type="field")
  # TODO: Need to replace, since using time-dep Neumann
  cellDomWhole.AssignBC()
  cellDomWhole.problem.d_eff = cellDomUnit.problem.d_eff
  print "Solving molecular unit cell using %s"% (meshFileCellular)
  solve_homogenized_whole(cellDomWhole.problem,cellDomUnit.problem,molDomUnit.problem,debug=debug)


# example 2.7 in Goel paper
##def GoelEx2p7():
##  ## micro
##  unitLength = np.array([0.1,0.1,1.0]) # length of microdomain in unit cell
##  diff = Constant(2.5) # isotropic
##
##  # our geometry needs to look something like the 'cross' in Fig 2
##  # creating this in Blender
##  # mesh is crap, doesn't work
##  meshFile = "goel.xml.gz"
##  problem = empty()
##  problem.mesh = Mesh(meshFile)
##
##  #mesh = UnitCube(6,6,6)
##  #problem.name ="test"
##  #problem.mesh = mesh
##
##  solve_homogeneous_unit(problem)
##  quit()
##
##
##  # note: authors take u0(x=[0.5,0.5,0,5]) = 0 to fix arb. const
##  # i think I can treat this like an initial condition and set some
##  # location in u0 to be zero.
##  u0 = Constant(0)  # this doesn't do anything, but could use as templat ealter
##  u_1 = interpolate(u0,V) # applies Expression u0 to FunctionSpace
##  # solve stuff
##  # u_1.assign(u)


##
## MAIN
##

if __name__ == "__main__":
  msg="script.py <arg>"
  remap = "none"

  #GoelEx2p7()

  import sys
  if len(sys.argv) < 1:
      raise RuntimeError(msg)


  # paths hardcoded insside
  #Debug()
  #Debug2()
  #quit()


  debug =0
  SolveHom(debug=debug)

