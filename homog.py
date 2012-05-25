#
# source config.bash
# from homog import *
# (cell,mol) = domost()

## TODO
# validate against COMSOL
# validate against Goel? (if mesh isn't too difficult)
# need to check on how to apply periodic BC 
# understand what BC to apply to problem 
# verify 'LCC' flux is being applied appopriately 
# Looks like I need to homogenize the fluxes as well (see Higgins) 
# add in smoluchowski component 

# Sub domain for Periodic boundary condition
# NOT PRESENTLY USED 
class PeriodicBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)

    def map(self, x, y):
        y[0] = x[0] - 1.0
        y[1] = x[1]


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

iso=1 


outbase= "/tmp/outs/"

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

# solve steady state diffusion problem with anisotropic diff constant
def build_steadystate(theDomain):
  problem = theDomain.problem

  # 3x3 diff matric 
  Dii  = Constant((problem.d_eff[0],problem.d_eff[1],problem.d_eff[2]))
  Dij = diag(Dii)  # for now, but could be anisotropic

  # build LHS
  u,v = TrialFunction(problem.V), TestFunction(problem.V)
  A = inner(Dij*grad(u), grad(v))*dx

  # init [not actually used in PDE soln]
  problem.x = Function(problem.V)
  problem.x.vector()[:] = parms.concInitial

  outname =outdir+problem.name+"_homogenized_stdy.pvd"
  print "Writing outputs to %s" % outname 
  out = File(outname) 

  # assign 
  problem.A = A
  problem.u = u
  problem.v = v
  problem.out = out

# 
# solves steady state diuffusion equation using flux, outerbounday conc
# provided by user 
def solve_steadystate(theDomain,outerconc,flux):
  problem = theDomain.problem 
  A = problem.A
  v = problem.v

  # source term indicative of flux betwen compartments
  n = FacetNormal(problem.mesh)
  L = flux*dot(n,problem.v)*ds
  
  boundaryConc = Constant((outerconc,outerconc,outerconc))
  theDomain.AssignBC( uBoundary=boundaryConc )

  x = Function(problem.V)
  solve(A==L,x,problem.bcs)

  problem.x.vector()[:] = x.vector()[:]

  # project and store
  #Vp = FunctionSpace(problem.mesh,"CG",1)
  #up = project(x[0], Vp)
  #problem.out << up
  problem.out << x 



# add in (time-dependent) flux condition
def buildRHSFlux(problem,t=0,j=0):
  print "Replace w real flux"
  mesh = problem.mesh
  n = FacetNormal(mesh)


  ## NOTE: I think it makes sense to keep this as simple as possible (e.g. no complicated fluxes or markers) 
  ## since mostly likely what I want is already in subcell. 
  flux = problem.dudn # hopefull the expression is getting updated with each time step 
  flux.t = t # update Expression 

  E = flux*dot(n,problem.v)*ds
  b  = assemble(E)

  # add in flux, if non-zero
  if(j!=0):
    E1 = j*dot(n,problem.v)*ds
    b1 = assemble(E)
    b += b1    

  return b

# builds matrices, etc, for time dependent solution 
def build_timedep(theDomain):
  problem = theDomain.problem 

  # 3x3 diff matric 
  print "Deff"
  print problem.d_eff
  #quit()
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
  b = buildRHSFlux(problem)
  # for multiple bcs (though currently we do not have any dirichlet)
  for bc in problem.bcs:
     bc.apply(A,b)


  # init cond
  u_n = Function(problem.V)
  xS = u_n.vector()
  xS[:] = parms.concInitial   
  problem.x = Function(problem.V) # not sure why I did this 
  problem.x.vector()[:] = xS[:]

  outname =outdir+problem.name+"_homogenized_tdep.pvd"
  print "Writing outputs to %s" % outname 
  out = File(outname) 

  # make assignments
  problem.K = K
  problem.A = A
  problem.M = M
  problem.b = b
  problem.xS = xS
  problem.out = out

# update matrices in time dep eqns
def update_timedep(theDomain,dt,flux):
  problem = theDomain.problem

  # LHS
  problem.A.assign(problem.K)                                              
  print "I don't think this term is correct for the 3x3 diff matrix"           
  problem.A *= parms.d*dt                                                    
  problem.A += problem.M         

  # RHS
  b1 = buildRHSFlux(problem,j=flux)
  problem.b=b1

  # add in mass matri
  problem.b += problem.M*problem.xS

  for bc in problem.bcs:
    bc.apply(problem.A,problem.b)

  ## solve and store 
  solve(problem.A,problem.xS,problem.b,"gmres","default")

  # store solution 
  problem.x.vector()[:] = problem.xS[:]

  # store results 
  #Vp = FunctionSpace(problem.mesh,"CG",1)
  #up = project(problem.x[0], Vp)
  #problem.out << up
  problem.out << problem.x

  
 
## Coupled problem (need to check on neumann confs)
# NOTE: assuming molecular domain is in steady state, so only solving
# cellular domain in time-dep fashion
def solve_homogenized_whole(wholecellDomain,unitcellDomain,unitmolDomain,type="field",debug=0):
  # get problems
  wholecell = wholecellDomain.problem
  unitcell = unitcellDomain.problem
  unitmol = unitmolDomain.problem


  ## molecular domain (steady state)
  print "Solving molecular steady state..."
  build_steadystate(unitmolDomain)

  ## Wholecell
  build_timedep(wholecellDomain)


  ## time parms 
  dt =parms.dt
  tStep = parms.tStep
  t = 0.
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
    k = 0.001
    hack = 0.5
    jcell = k*(wholecell.conc - hack*unitmol.conc)
    print "wholecell: %f" % wholecell.conc
    print "unitmol: %f" % unitmol.conc
    jmol  = -jcell

    # TODO add scale factors (are these consistent)
    jcell*= unitcell.surfaceArea/ unitcell.gamma
    jmol *= unitmol.surfaceArea/ unitmol.gamma
    print "jcell %f" % jcell
    
    # solve cell
    # JOHAN
    print "t %f" % t
    t  += float(dt)
    
    print "Need to add in coupling between molecular and cellular"
    update_timedep(wholecellDomain,dt,jcell)

    #  TODOGoel: how/where is the surface defined here? It seems like the boundaries within
    # the unit cell are 'embedded' into the large cell description, so are no longer
    # boundaries
    print" TODO - note: ds needs to be limited to boundary between cell and microdomains [eq use ds(2)]"
    solve_steadystate(unitmolDomain,wholecell.conc,jmol)


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


##
## Example calls (should go in its own file)
##


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

  #molDomUnit.problem.d_eff = np.array([1.,1.,1.])
  if(iso==1):
    print "WARNING: I am CHEATING by overriding anistropic diff const"
    dcheatmol= np.array([1.,1.,1.])
    dcheatcell = np.array([1.,1.,1.])
  elif(iso==2):
    print "WARNING: I am CHEATING by overriding anistropic diff const"
    dcheatmol= np.array([0.1,0.1,10.])
    dcheatcell = np.array([10.,0.1,0.1])

  molDomUnit.problem.d_eff = dcheatmol
  cellDomUnit.problem.d_eff = dcheatcell
  print dcheatmol
  print dcheatcell

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

def ComsolEx():

  molDomUnit = DefaultUnitDomain(type="field")
  molDomUnit.Setup()
  molDomUnit.AssignBC()

  #mesh1= UnitSphere(5,5,5)
  # make same size as comsol ex. 
  molDomUnit.problem.mesh.coordinates()[:]*0.5005 

  # not sure if cube is even used in the model - looks like just a sphere here
  #mesh2 = UnitCube(5,5,5)
  # Assumed in params()

  #convective flux source of -1?? - ah, due to (del^2 u + I)

  # NOTE: need to modify to handle PBC for field, since constraint loation is differnet for each solution vector component
  print "Need to enforce periodic bc for field" 
  # looks like there's some kind of PBC at either side of the sphere (continuity)
  # define 'periodic boundary' at R > 0.5 (hence using a radius of 0.5005)
  pbc = PeriodicBoundary()
  bc1 = PeriodicBC(V, pbc)
  print "need to impose PBC for FLUXes as well"
#
  identifies pt at 0,0,0 
  # zero flux applied to all sphere boundaries 

  ## soln 

  #???CalcFractionalVolumes(cellDomUnit,molDomUnit)

  solve_homogeneous_unit(cellDomUnit,type="field")
  


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
  msg="script.py <name>"
  remap = "none"

  #GoelEx2p7()

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)
  
  outdir = outbase+'/'+sys.argv[1]+'/'
  print "Writing outputs to %s" % outdir

  # paths hardcoded insside
  ComsolEx()
  #Debug2()
  #quit()
  sys.argv[1]
  if(sys.argv[1]=="isocheat"):
    iso=1
  elif(sys.argv[1]=="nonisocheat"):
    iso=2


  # globular case
  #SolveGlobularHomog(debug=debug)

  # TnC/cylindrical case
  #SolveMyofilamentHomog()

