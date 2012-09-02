#
# source config.bash
# from homog import *
# (cell,mol) = domost()

## TODO
# DONE add in PBC for cell
# use unit cell geometry that is based on the molec geom (right now its sphereically symm)  
# validate against COMSOL
# validate against Goel? (if mesh isn't too difficult)
# verify 'LCC' flux is being applied appopriately 
# Looks like I need to homogenize the fluxes as well (see Higgins) 



from dolfin import *
import numpy as np
from params import *
import matplotlib.pyplot as plt

# classes
from CellularDomain import *
from DefaultUnitDomain import *
from CellularUnitDomain import *
from MolecularUnitDomain import *

from CellularUnitDomain_TnC import *

import scalar
import field

##
## VAR
## 
debug=0
smolMode = "false" # tells program that we want to solve the smol equation for the molec domain
case ="none"
molGamer = 0 
validationMode = 0
#outbase= "/tmp/outs/"
outbase="/home/huskeypm/scratch/homogouts/"

# rocce
cellPrefix="none"
wholeCellPrefix="none"
molPrefix="none"
root = "/home/huskeypm/scratch/homog/"

# vm
#root = "/home/huskeypm/localTemp/"
 

class empty:pass

## solv. homog cell
def solve_homogeneous_unit(domain,type="field",debug=0,smolMode="false"):
  problem = domain.problem

  ## debug mode
  if(debug==1):
  #if(1):
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
    scalar.solveHomog(domain)
    D_eff = scalar.compute_eff_diff(domain)
    d_eff = np.array([D_eff,D_eff,D_eff])

  # using vector fields 
  elif (type=="field"):
    field.solveHomog(domain,smolMode=smolMode)
    d_eff = field.compute_eff_diff(domain)

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

  outname =outbase+problem.name+"_homogenized_stdy.pvd"
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
# Specific to molecule right now 
def solve_molecular_steadystate(molDomain,outerconc,flux):
  problem = molDomain.problem 
  A = problem.A
  v = problem.v

  # source term indicative of flux betwen compartments
  n = FacetNormal(problem.mesh)
  L = flux*dot(n,problem.v)*ds(molDomain.markerOuterBoundary)
  
  boundaryConc = Constant((outerconc,outerconc,outerconc))
  molDomain.AssignBC( uBoundary=boundaryConc )

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
def build_timedep(theDomain,tag=""):
  problem = theDomain.problem 

  # 3x3 diff matric 
  print "Using Deff for wholecell"
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

  outname =outbase+tag+problem.name+"_homogenized_tdep.pvd"
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
def solve_homogenized_whole(wholecellDomain,unitcellDomain,unitmolDomain,type="field",tag="",debug=0,wholeCellOnly=0):
  # get problems
  wholecell = wholecellDomain.problem
  unitcell = unitcellDomain.problem
  unitmol = unitmolDomain.problem


  ## molecular domain (steady state)
  if(wholeCellOnly!=0):
    print "Solving molecular steady state..."
    build_steadystate(unitmolDomain)

  ## Wholecell
  build_timedep(wholecellDomain,tag=tag)
  wholecell.conc = parms.concInitial


  ## time parms 
  dt =parms.dt
  tStep = parms.tStep
  t = 0.
  tMax = tStep * dt


  ## iterate
  while (t < tMax):

    ## TODO check that flux is correct
    # assume simple flux between compartments
    # CHeck validation routine if adding this back in (cellOnly case will fail here) 
    if(wholeCellOnly!=0):
      field.CalcConc(unitmolDomain)
    else:
      unitmol.conc=0.001

    # TODO - flux between unitmol and wholecell domain (where is the surface here, since this
    # concept seems tobe relevant only for unit cell??) 
    # NOTE: I'm not so sure I need to have a flux between molecular domain and cellular domain, since
    # we assume that molecular domain is at steady state and Dirichlet on boundary is equal to cellular conc
    k = 100
    hack = 0.5
    jcell = k*(wholecell.conc - hack*unitmol.conc)
    print "wholecell: %f" % wholecell.conc
    print "unitmol: %f" % unitmol.conc
    jmol  = -jcell

    # TODO add scale factors (are these consistent)
    if(wholeCellOnly!=0):
      jcell*= unitcell.surfaceArea/ unitcell.gamma
      jmol *= unitmol.surfaceArea/ unitmol.gamma

    # some other flux
    flux = 100 * np.exp(-t/10)
    jcell += flux


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
    print "WARNING: skipping molec stdy st until boundary marking issue fixed"
    #DEBUGsolve_molecular_steadystate(unitmolDomain,wholecell.conc,jmol)


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


# Solves homogenized equations for globular protein embedded inside spherical cellular subdomain 
def SolveHomogSystem(debug=0,\
  root="./",\
#  cellPrefix="cell",molPrefix="mol",wholeCellPrefix="multi_clustered",\
  cellPrefix="none",molPrefix="none",wholeCellPrefix="none",\
# use effective diffusion constant from molecular domain 
  useMoldeff=1,\
  smolMode = "false",\
# is molecule from Gamer?
  molGamer=1,\
# for passing other options
  option="",\
  tag=""):   # optional tag for file naming 
   

  ## Setup
  #root = "/home/huskeypm/scratch/homog/mol/"

  results = empty()

  # celular domain
  #if(debug==0):
  if(cellPrefix!="none"):
    meshFileOuter = root+cellPrefix+"_mesh.xml.gz"
    subdomainFileOuter = root+cellPrefix+"_subdomains.xml.gz"
    cellDomUnit = CellularUnitDomain(meshFileOuter,subdomainFileOuter,\
      type="field",outpath=outbase,name=tag+cellPrefix)
    cellDomUnit.Setup()
    cellDomUnit.AssignBC()
    results.cellDomUnit = cellDomUnit
  else:
    cellDomUnit = empty()
    cellDomUnit.problem = empty()


  # molecular domain
  if(molPrefix!="none"): 
    meshFileInner = root+molPrefix+"_mesh.xml.gz"
    subdomainFileInner = root+molPrefix+"_subdomains.xml.gz"
    potentialFileInner = root+molPrefix+"_values.xml.gz"
    molDomUnit = MolecularUnitDomain(meshFileInner,subdomainFileInner,\
      filePotential = potentialFileInner,type="field",gamer=molGamer,\
      outpath=outbase,name=tag+molPrefix)
    molDomUnit.problem.smolMode = smolMode
    # 
    if(option=="troponin"):
      molDomUnit.problem.scale = np.array([0.9*np.sqrt(2)/2.,0.9*np.sqrt(2)/2.,0.9])
      print "WARNING: I AM USING A SCALE FACTOR HERE TO ASSIGN BOUNDARY "
      print molDomUnit.problem.scale

    molDomUnit.Setup()
    molDomUnit.AssignBC()
    results.molDomUnit = molDomUnit
  else:
    molDomUnit = empty()
    molDomUnit.problem = empty()


  # get fractional volume 
  #if(debug==0):
  if(molPrefix!="none" and cellPrefix!="none"): 
    CalcFractionalVolumes(cellDomUnit,molDomUnit)

  ## Solve unit cell problems  
  #if(debug==0):
  if(cellPrefix!="none"): 
    print "Solving cellular unit cell using %s" % meshFileOuter
    solve_homogeneous_unit(cellDomUnit,type="field",debug=debug)
  else:
    print "WARNING: Using stored d_eff for cell" 
    cellDomUnit.problem.d_eff = np.array([1,1,1])

  if(molPrefix!="none"): 
    print "Solving molecular unit cell using %s"% meshFileInner
    molDomUnit.smolMode = smolMode
    solve_homogeneous_unit(molDomUnit,type="field",debug=debug,smolMode=smolMode)
  else:
    print "WARNING: Using stored d_eff for mol"   
    molDomUnit.problem.d_eff = np.array([1,1,1])

     


  ## whole cell solutions
  # solve on actual cell domain
  if(wholeCellPrefix!="none"):
    meshFileCellular= root+wholeCellPrefix+"_mesh.xml.gz"
    subdomainFileCellular= root+wholeCellPrefix+"_subdomains.xml.gz"
    if(debug==0):
      cellDomWhole = CellularDomain(meshFileCellular,subdomainFileCellular,type="field")
    else:
      cellDomWhole = DefaultUnitDomain(type="field")
  
    cellDomWhole.Setup()
    cellDomWhole.AssignBC()
    results.cellDomWhole= cellDomWhole
  
  
    # assign diff const 
    if(useMoldeff==1):
      print "WARNING: For now will assume cellular diffusion dominated by molec dom"
      cellDomWhole.problem.d_eff = molDomUnit.problem.d_eff
    else:
      cellDomWhole.problem.d_eff = cellDomUnit.problem.d_eff
  
    # solve 
    if(debug==0):
      print "Solving macro equations using %s"% (meshFileCellular)
      solve_homogenized_whole(cellDomWhole,cellDomUnit,molDomUnit,tag=tag,debug=debug)


  ## return 
  return results

# ValidationSphere case for simple charged sphere 
def ValidationSphere():
  root = "/home/huskeypm/scratch/validation/sphere/"
  molPrefix = "sphere"
  #molDomUnit.problem.mesh.coordinates()[:]*0.5005 

  ## simple sphere 
  molDomUnit = DefaultUnitDomain()
  molDomUnit.Setup()
  molDomUnit.AssignBC()
  problem = molDomUnit.problem
  problem.pmf = Function( FunctionSpace(problem.mesh,"CG",1))
  problem.pmf.vector()[:] = np.arange(729)/729. - 0.5
  #problem.pmf.vector()[:] = 0
  File("test.pvd") << problem.pmf
  #problem.pmf.vector()[:] = 0
  smolMode = "true"
  smolMode = "false"
  solve_homogeneous_unit(molDomUnit,type="field",debug=debug,smolMode=smolMode)
  quit()

  ## gamer sphere 
  print "No electro" 
  smolMode = "false" # tells program that we want to solve the smol equation for the molec domain
  noelectroResults = SolveHomogSystem(debug=debug,\
    root=root,\
    molPrefix=molPrefix,
    smolMode = smolMode,\
    molGamer=0)

  print "With electro"
  smolMode = "true" # tells program that we want to solve the smol equation for the molec domain
  electroResults = SolveHomogSystem(debug=debug,\
    root=root,\
    molPrefix=molPrefix,
    smolMode = smolMode,\
    molGamer=0)

  # compare diff const
  d_eff_noelectro = noelectroResults.molDomUnit.problem.d_eff
  nd_eff_noelectro = d_eff_noelectro / np.linalg.norm(d_eff_noelectro)
  d_eff_electro = electroResults.molDomUnit.problem.d_eff
  nd_eff_electro = d_eff_electro / np.linalg.norm(d_eff_electro)
  print "Deff (No Electro) ", nd_eff_noelectro
  print "Deff (Electro) ", nd_eff_electro

def ValidationLattice():
    allsummary=[]
    mode="valid"
    #python ~/localTemp/srcs/homogenization/homog.py -case custom  -molPrefix 1p50 -molGamer
    molPrefix = "1p50"
    #molPrefix = "test"
    molPrefixes = ["0p50","0p75","1p01","1p25","1p50"]
    molDists    = 2 - np.array([0.50,0.75,1.01,1.25,1.50])
    norms = np.zeros( len(molDists) )

    
    molGamer = 1
    for i,molPrefix in enumerate(molPrefixes):
      print "molPrefix"
      results = SolveHomogSystem(debug=debug,\
        root="./validation/",\
        cellPrefix="none", molPrefix=molPrefix,wholeCellPrefix="none",\
        smolMode = "false",\
        molGamer=molGamer,
        tag=mode)

      norms[i]  = np.linalg.norm(results.molDomUnit.problem.d_eff)
      r=results.molDomUnit.problem.d_eff    
      r = r/np.linalg.norm(r)
      summary = "%s & Deff (%e,%e,%e) \\ \n" % (mode,r[0],r[1],r[2])
      allsummary.append(summary)

    # plot 
    plt.figure()
    plt.plot(molDists,norms)
    plt.scatter(molDists,norms)
    plt.ylabel("|D| [$m^2/s$]")
    plt.xlabel("Edge length [m]")
    plt.title("Diffusivity versus molecule size")
    f=plt.gcf()
    f.savefig("Boxes.png")

    return (molDists,norms)


# Paper validation, fig gen
def ValidationPaper(mode="all"):
  allsummary=[]

  if(mode == "lattice" or mode == "all"):
    ValidationLattice()

  if(mode=="troponinNoChg" or mode =="all"):
    molPrefix = "troponin"
    results = SolveHomogSystem(debug=debug,\
      root="/home/huskeypm/scratch/homog/mol/",\
      cellPrefix="none", molPrefix=molPrefix,wholeCellPrefix="none",\
      smolMode = "false",
      molGamer=0,option="troponin",\
      tag=mode)

    r=results.molDomUnit.problem.d_eff    
    r = r/np.linalg.norm(r)
    summary = "%s & Deff (%e,%e,%e) \\ \n" % (mode,r[0],r[1],r[2])
    allsummary.append(summary)

  if(mode=="troponinWChg" or mode =="all"):
    molPrefix = "troponin"
    results = SolveHomogSystem(debug=debug,\
      root="/home/huskeypm/scratch/homog/mol/",\
      cellPrefix="none", molPrefix=molPrefix,wholeCellPrefix="none",\
      smolMode = "true",
      molGamer=0,option="troponin",\
      tag=mode)

    r=results.molDomUnit.problem.d_eff    
    r = r/np.linalg.norm(r)
    summary = "%s & Deff (%e,%e,%e) \\ \n" % (mode,r[0],r[1],r[2])
    allsummary.append(summary)

  if(mode=="cellOnly"):
    #parms.tStep = 2
    wholeCellPrefix="multi_clustered"
    results = SolveHomogSystem(debug=debug,\
      root="/home/huskeypm/scratch/homog/mol/",\
      cellPrefix="none", molPrefix="none",wholeCellPrefix=wholeCellPrefix,\
      smolMode = "false",
      molGamer=0,\
      tag=mode)

      # no need to print diff 

  if(mode=="totalNoChg" or mode=="all"): 
    #parms.tStep = 2
    cellPrefix="cell"
    wholeCellPrefix="multi_clustered"
    molPrefix = "troponin"
    results = SolveHomogSystem(debug=debug,\
      root="/home/huskeypm/scratch/homog/mol/",\
      cellPrefix=cellPrefix, molPrefix=molPrefix,wholeCellPrefix=wholeCellPrefix,\
      smolMode = "false",
      molGamer=0,option="troponin",\
      tag=mode)

    r=results.molDomUnit.problem.d_eff    
    r = r/np.linalg.norm(r)
    summary = "%s & Deff (%e,%e,%e) \\ \n" % (mode,r[0],r[1],r[2])
    allsummary.append(summary)

  if(mode=="totalWChg" or mode =="all"):
    cellPrefix="cell"
    wholeCellPrefix="multi_clustered"
    molPrefix = "troponin"
    results = SolveHomogSystem(debug=debug,\
      root="/home/huskeypm/scratch/homog/mol/",\
      cellPrefix=cellPrefix, molPrefix=molPrefix,wholeCellPrefix=wholeCellPrefix,\
      smolMode = "true",
      molGamer=0,option="troponin",\
      tag=mode)

    r=results.molDomUnit.problem.d_eff    
    r = r/np.linalg.norm(r)
    summary = "%s & Deff (%e,%e,%e) \\ \n" % (mode,r[0],r[1],r[2])
    allsummary.append(summary)


  print allsummary





  


##
## MAIN
##



if __name__ == "__main__":
  msg="""
\nPurpose: 
  run homogenized problem 
 
Usage:
  homog.py -case myofilament/globular/validationSphere/custom <-smol> <-molGamer>
           <-molPrefix molPrefix> <-validation all/etc>

  where 
    -smol - run molecular domain with electrostatics
    -molGamer - molecule was prepared with Gamer

Notes:
"""

  #GoelEx2p7()

  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)
  
  print "Writing outputs to %s" % outbase
  #print "WARNING: files arent writing!"

  for i,arg in enumerate(sys.argv):
    if(arg=="-smol"):
      print "Ensabling electrostatics contribution" 
      smolMode = "true"

    if(arg=="-case"):
      case = sys.argv[i+1]

    if(arg=="-molPrefix"):
      molPrefix = sys.argv[i+1]

    if(arg=="-molGamer"):
      molGamer = 1

    if(arg=="-validation"):
      validationMode=sys.argv[i+1]


  #
  # Validation 
  #  
  if(validationMode!=0):          
    ValidationPaper(mode=validationMode)
    print "Run paraview.py to view results"
    quit()


  #
  # Normal 
  #  


  #Debug2()
  #quit()

  # ovoerride
  #debug =1
  #cellPrefix = ""
  #molPrefix = "molecular_volmesh"
  #root = "/home/huskeypm/bin/grids/"
  
  #case = "globular"
  #case = "myofilament"
  #case = "validation"


  # globular case
  if(case=="globular"):
    #cellPrefix="mol/cell"
    #wholeCellPrefix="mol/multi_clustered"
    #molPrefix="120529_homog/1CID"
    molPrefix = "mol/1CID"
    SolveHomogSystem(debug=debug,\
      root=root,\
      cellPrefix=cellPrefix, molPrefix=molPrefix,wholeCellPrefix=wholeCellPrefix,\
      molGamer=1)
  
  # TnC/cylindrical case
  elif(case=="myofilament"):
    cellPrefix="mol/cell"
    wholeCellPrefix="mol/multi_clustered"
    molPrefix = "mol/troponin"
    SolveHomogSystem(debug=debug,\
      root=root,\
      cellPrefix=cellPrefix, molPrefix=molPrefix,wholeCellPrefix=wholeCellPrefix,\
      molGamer=molGamer,option="troponin")


  elif(case=="validationSphere"):
    ValidationSphere()

  elif(case=="custom"):
    root = "./"
    SolveHomogSystem(debug=debug,\
      root=root,\
      cellPrefix=cellPrefix, molPrefix=molPrefix,wholeCellPrefix=wholeCellPrefix,\
      smolMode = smolMode,
      molGamer=molGamer)

  else:
    msg = "Case " + case + " not understood"   
    raise RuntimeError(msg)


