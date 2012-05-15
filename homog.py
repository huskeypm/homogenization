#
# source config.bash
# from homog import *
# (cell,mol) = domost()

## TODO
# validate against COMSOL
# validate against Goel? (if mesh isn't too difficult)
# understand what BC to apply to problem 
# verify 'LCC' flux is being applied appopriately 
# Looks like I need to homogenize the fluxes as well (see Higgins) 
# add in smoluchowski component 

from dolfin import *
import numpy as np
from params import *

# classes
from CellularDomain import *
from DefaultUnitDomain import *
from CellularUnitDomain import *
from MolecularUnitDomain import *

from CellularUnitDomain_TnC import *

import scalar
import field



class empty:pass


## solv. homog cell
def solve_homogeneous_unit(domain,type="field",debug=0):
  problem = domain.problem

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

# builds matrices, etc, for time dependent solution 
def build_timedep(theDomain):
  problem = theDomain.problem 

  # 3x3 diff matric 
  Dii  = Constant((problem.d_eff[0],problem.d_eff[1],problem.d_eff[2]))
  Dij = diag(Dii)  # for now, but could be anisotropic

  # Assembly of the K, M and A matrices
  u,v = TrialFunction(problem.V), TestFunction(problem.V)
  problem.u = u
  problem.v = v
  a =inner(Dij * grad(u), grad(v))*dx
  K = assemble(a,mesh=problem.mesh)
  A = K.copy()
  # TODO check on this
  M = assemble(inner(u,v)*dx,mesh=problem.mesh)

  # TODO check on this 
  b = getb(problem)
  # for multiple bcs (though currently we do not have any dirichlet)
  for bc in problem.bcs:
     bc.apply(A,b)


  # init cond
  u_n = Function(problem.V)
  x = u_n.vector()
  x[:] = 0.1 # initial conc  
  problem.x = Function(problem.V) # not sure why I did this 
  problem.x.vector()[:] = x[:]

  out  = File(problem.name+"_homogenized.pvd")

  # make assignments
  problem.K = K
  problem.A = A
  problem.M = M
  problem.b = b
  problem.out
  


 
## Coupled problem (need to check on neumann confs)
# NOTE: assuming molecular domain is in steady state, so only solving
# cellular domain in time-dep fashion
def solve_homogenized_whole(wholecellDomain,unitcellDomain,unitmolDomain,type="field",debug=0):
  # get problems
  wholecell = wholecellDomain.problem
  unitcell = unitcellDomain.problem
  unitmol = unitmolDomain.problem


  # store output 

  ## molecular domain (steady state)
  # make Diff matrix
  print "Solving molecular steady state..."
  #Dii  = Constant((wholecell.d_eff[0],wholecell.d_eff[1],wholecell.d_eff[2]))
  #Dij = diag(Dii)  # for now, but could be anisotropic
  #u,v = TrialFunction(unitmol.V), TestFunction(unitmol.V)
  #A = inner(Dij*grad(u), grad(v))*dx
  #L = (beta  * v)*dx
  #x = Function(unitmol.V)
  d = 1
  Dii  = Constant((d,d,d))
  Dij = diag(Dii)  # for now, but could be anisotropic
  
  u,v = TrialFunction(unitmol.V), TestFunction(unitmol.V)
  A = inner(Dij*grad(u), grad(v))*dx
  # source term indicative of flux betwen compartments
  c = Constant((1,1,1))
  L = inner(c,v)*dx
  
  print "Should I use avg cell conc as dirichlet?"
  print "specific to moleculare case"
  unitmolDomain.AssignBC( uBoundary=Constant((1.,1.,1.)) )
  #bc0 = DirichletBC(unitmol.V, Constant((0.,0.,0.)), unitmol.subdomains,1)
  #bc1 = DirichletBC(unitmol.V, Constant((1.,1.,1.)), unitmol.subdomains,5)
  #bcs = [bc0,bc1]
  x = Function(unitmol.V)
  solve(A==L,x,unitmol.bcs)
  File(unitmol.name+"_homogenized.pvd") << x
  print "TODO: start w initial conc and add source term. Not sure what BC"


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
    # NOTE: I'm not so sure I need to have a flux between molecular domain and cellular domain, since
    # we assume that molecular domain is at steady state and Dirichlet on boundary is equal to cellular conc
    k = 1
    hack = 0.5
    jcell = k*(wholecell.conc - hack*unitmol.conc)
    print "wholecell: %f" % wholecell.conc
    print "unitmol: %f" % unitmol.conc
    jmol  = -jcell

    # TODO add scale factors (are these consistent)
    #print "WARNING: cmtd out"
    #jcell*= unitcell.beta / unitcell.gamma
    #jmol *= unitmol.beta / unitmol.gamma
    
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
    print "probably either want to project or plot single component"
    out << wholecell.x


    # TODO? solv unitmol
    #F = unitmol.diff_eff * grad(unitmol.u) * grad(unitmol.v)
    #F += unitmol.diff_eff( junitmol * unitmol.v)
    #M = 0
    #solve(F==M,unitmol.u)
    #write



  print "Finished!"

def CalcFractionalVolumes(cellDomUnit,molDomUnit):
  cellProblem = cellDomUnit.problem
  molProblem = molDomUnit.problem
  volUnitCell = cellProblem.volume + molProblem.volume 
  cellProblem.gamma = cellProblem.volume/volUnitCell
  cellProblem.volUnitCell = volUnitCell
  print "cell frac vol: %f" % cellProblem.gamma
  molProblem.gamma = molProblem.volume/volUnitCell
  molProblem.volUnitCell = volUnitCell
  print "mol frac vol: %f" % molProblem.gamma



def Debug():
  cellDomUnit = DefaultUnitDomain(type="field")
  cellDomUnit.Setup()
  cellDomUnit.AssignBC()

  solve_homogeneous_unit(cellDomUnit,type="field")

# This function should evolve into a general protocol for solvign the micro/macro equations 
def Debug2():
  ## microdomain problems
  cellDomUnit = DefaultUnitDomain(type="field")
  cellDomUnit.Setup()
  cellDomUnit.AssignBC()


  molDomUnit = DefaultUnitDomain(type="field")
  molDomUnit.Setup()
  molDomUnit.AssignBC()

  CalcFractionalVolumes(cellDomUnit,molDomUnit)

  solve_homogeneous_unit(cellDomUnit,type="field")
  solve_homogeneous_unit(molDomUnit,type="field")


  ## macrodomain problems 
  cellDomWhole = DefaultUnitDomain(type="field")
  cellDomWhole.Setup()
  cellDomWhole.AssignBC()
  cellDomWhole.problem.d_eff = cellDomUnit.problem.d_eff
  solve_homogenized_whole(cellDomWhole,cellDomUnit,molDomUnit)




# Solve homogenized equations for myofilament embedded in rectilinear cellular domain
# See notes from 120514_grids.tex as what we are doing here is a bit counterintuitive 
def SolveMyofilamentHomog(root="./",debug=0):

  ## Setup
  root = "/home/huskeypm/scratch/homog/mol/"

  # celular domain
  prefix = "cell"
  cellDomUnit = CellularUnitDomain_TnC(type="field")
  cellDomUnit.Setup()
  cellDomUnit.AssignBC()

  # molecular domain
  prefix = "troponin"
  meshFileInner = root+prefix+"_mesh.xml.gz"
  subdomainFileInner = root+prefix+"_subdomains.xml.gz"
  molDomUnit = MolecularUnitDomain(meshFileInner,subdomainFileInner,type="field")
  molDomUnit.Setup()
  molDomUnit.AssignBC()

  # get fractional volume 
  CalcFractionalVolumes(cellDomUnit,molDomUnit)

  ## Solve unit cell problems  
  print "Solving cellular unit cell "
  solve_homogeneous_unit(cellDomUnit,type="field",debug=debug)
  print "Solving molecular unit cell using %s"% meshFileInner
  solve_homogeneous_unit(molDomUnit,type="field",debug=debug)

  print "WARNING: I am CHEATING by overriding anistropic diff const"
  molDomUnit.problem.d_eff = np.array([1.,0.1,0.1])
  cellDomUnit.problem.d_eff = np.array([1.,0.1,0.1])

  ## whole cell solutions
  # solve on actual cell domain
  prefix = "multi_clustered"
  meshFileCellular= root+prefix+"_mesh.xml.gz"
  subdomainFileCellular= root+prefix+"_subdomains.xml.gz"
  if(debug==0):
    cellDomWhole = CellularDomain(meshFileCellular,subdomainFileCellular,type="field")
  else:
    cellDomWhole = DefaultUnitDomain(type="field")

  cellDomWhole.Setup()
  cellDomWhole.AssignBC()
  cellDomWhole.problem.d_eff = cellDomUnit.problem.d_eff
  print "Solving molecular unit cell using %s"% (meshFileCellular)
  solve_homogenized_whole(cellDomWhole,cellDomUnit,molDomUnit,debug=debug)

# Solves homogenized equations for globular protein embedded inside spherical cellular subdomain 
def SolveGlobularHomog(root="./",debug=0):

  ## Setup
  root = "/home/huskeypm/scratch/homog/mol/"

  # celular domain
  prefix = "cell"
  meshFileOuter = root+prefix+"_mesh.xml.gz"
  subdomainFileOuter = root+prefix+"_subdomains.xml.gz"
  cellDomUnit = CellularUnitDomain(meshFileOuter,subdomainFileOuter,type="field")
  cellDomUnit.Setup()
  cellDomUnit.AssignBC()

  # molecular domain
  prefix = "mol"
  meshFileInner = root+prefix+"_mesh.xml.gz"
  subdomainFileInner = root+prefix+"_subdomains.xml.gz"
  molDomUnit = MolecularUnitDomain(meshFileInner,subdomainFileInner,type="field")
  molDomUnit.Setup()
  molDomUnit.AssignBC()

  # get fractional volume 
  CalcFractionalVolumes(cellDomUnit,molDomUnit)


  ## Solve unit cell problems  
  print "Solving cellular unit cell using %s" % meshFileOuter
  solve_homogeneous_unit(cellDomUnit,type="field",debug=debug)
  print "Solving molecular unit cell using %s"% meshFileInner
  solve_homogeneous_unit(molDomUnit,type="field",debug=debug)


  ## whole cell solutions
  # solve on actual cell domain
  prefix = "multi_clustered"
  meshFileCellular= root+prefix+"_mesh.xml.gz"
  subdomainFileCellular= root+prefix+"_subdomains.xml.gz"
  if(debug==0):
    cellDomWhole = CellularDomain(meshFileCellular,subdomainFileCellular,type="field")
  else:
    cellDomWhole = DefaultUnitDomain(type="field")

  cellDomWhole.Setup()
  cellDomWhole.AssignBC()
  cellDomWhole.problem.d_eff = cellDomUnit.problem.d_eff
  print "Solving macro equations using %s"% (meshFileCellular)
  solve_homogenized_whole(cellDomWhole,cellDomUnit,molDomUnit,debug=debug)


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
  Debug2()
  quit()


  debug =0
  # globular case
  #SolveGlobularHomog(debug=debug)

  # TnC/cylindrical case
  SolveMyofilamentHomog()

