# Defines boundary conditions, etc for microdomain
from dolfin import *
from params import *
from Domain import *
from util import * 
import smol

markerActiveSite = 1
markerMolecularBoundary =4
markerOuterBoundary=5
q = 2.0 # charge Ca2+ 

EPS = 1.e-10

#PKH2class LeftRightBoundary(SubDomain):
#PKH2        def inside(self, x, on_boundary):
#PKH2            # find v1x
#PKH2            problem = self.problem 
#PKH2            return tuple(x) in problem.targetsx
#PKH2    
#PKH2    
#PKH2        # field component 1
#PKH2        def map(self, x, y):
#PKH2            problem = self.problem 
#PKH2            y[:] = problem.vert_mapx.get(tuple(x), x)
#PKH2    
#PKH2class BackFrontDomain(SubDomain):
#PKH2        def inside(self, x, on_boundary):
#PKH2            # find v1x
#PKH2            problem = self.problem 
#PKH2            return tuple(x) in problem.targetsy
#PKH2    
#PKH2    
#PKH2        # field component 1
#PKH2        def map(self, x, y):
#PKH2            problem = self.problem 
#PKH2            y[:] = problem.vert_mapy.get(tuple(x), x)
#PKH2    
#PKH2    
#PKH2class TopBottomDomain(SubDomain):
#PKH2        def inside(self, x, on_boundary):
#PKH2            # find v1x
#PKH2            problem = self.problem 
#PKH2            return tuple(x) in problem.targetsz
#PKH2    
#PKH2    
#PKH2        # field component 1
#PKH2        def map(self, x, y):
#PKH2            problem = self.problem 
#PKH2            y[:] = problem.vert_mapz.get(tuple(x), x)
#PKH2    
#PKH2    # Sub domain for Dirichlet boundary condition
#PKH2class CenterDomain(SubDomain):
#PKH2        def inside(self, x, in_boundary):
#PKH2            problem = self.problem 
#PKH2            return all(near(x[i], problem.center_coord[i], EPS) for i in range(problem.nDims))


class MolecularUnitDomain(Domain):
  def __init__(self,fileMesh,fileSubdomains,\
    # scalar (wrong) or field
    type,\
    # doe mesh come from gamer?
    gamer=1\
    ):
    super(MolecularUnitDomain,self).__init__(type)

    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.init = 1
    #print "Enforcing gamer==1"
    self.gamer = gamer
    problem.name = "Molecular"
    self.markerOuterBoundary = markerOuterBoundary

  def Setup(self):
    # mesh
    problem = self.problem
    mesh = Mesh(problem.fileMesh)
    problem.mesh = mesh
    # mesh is in A, so converting to um
    mesh.coordinates()[:] *= parms.ANG_TO_UM

    utilObj = util(problem)
    utilObj.GeometryInitializations(mesh)
    utilObj.DefinePBCMappings(mesh)

    # rotate to align z axis of molecule with x axis of myocyte
    print "WARNING: need to rotate mesh to align with myocyte. (my code seemed to cause PETSC error)"
    #x = mesh.coordinates()[:,0]
    #z = mesh.coordinates()[:,2]
    #mesh.coordinates()[:,0] = -1 * z[:]
    #mesh.coordinates()[:,2] = x[:]
    #V = FunctionSpace(mesh,"CG",1)
    #v = Function(V)
    #File("test.pvd") << v
    #quit()


    if(self.type=="scalar"):
        problem.V = FunctionSpace(mesh,"CG",1)
    elif(self.type=="field"):
        problem.V = VectorFunctionSpace(mesh,"CG",1)

    problem.subdomains = MeshFunction(
      "uint", mesh, problem.fileSubdomains)
    #self.markers = [1,4,5]

    # load ESP 
    if(0):
      print "PUT IN SEPARATE FUNCTION"
      Vtemp = VectorFunctionSpace(mesh,"CG",1)
      psi = Function(Vtemp,problem.filePotential)
      pmf = smol.ElectrostaticPMF(problem,psi,q=q, V=Vtemp)
      problem.pmf = Function(mesh.V)
      # assign to each cpomponent 
      print "WARNING: I am not sure if ESP is applied correctly here"
      problem.pmf.vector()[0,:] = pmf
      problem.pmf.vector()[1,:] = pmf
      problem.pmf.vector()[2,:] = pmf
    

    # geom
    self.CalcGeom(problem)


  # bcs
  def AssignBC(self,uBoundary=0):
    problem = self.problem

    #print "Probably don't have the right BCs yet"
    # TODO: might need to handle BC at active site as a mixed boundary
    # condition to take into account saturation
    if(self.type=="scalar"):
        print "REPLACE THIS AS DEBUG OPTION. scalar approach is nonsensical"
        quit()
        u0 = Constant(0.)
        u1 = Constant(1.)
    elif(self.type=="field"):
        u0 = Constant((0.,0,0.))
        u1 = Constant((1.,1.,1.))

    # use user-provided BC instead  
    if(uBoundary != 0):
      u1 = uBoundary

    ## PKH WAS 
    ## molec boundary
    #bc0 = DirichletBC(problem.V,u0,problem.subdomains,markerActiveSite)
    ## outer boundary
    #bc1 = DirichletBC(problem.V,u1,problem.subdomains,markerOuterBoundary)
    #problem.bcs = [bc0,bc1]

  
    ## Nested classes for handling periodic BCs
    # TODO: would like to embed this in Domain.py, not sure how to do this
#    class LeftRightBoundary(SubDomain):
#        def inside(self, x, on_boundary):
#            # find v1x
#            return tuple(x) in problem.targetsx
#    
#    
#        # field component 1
#        def map(self, x, y):
#            y[:] = problem.vert_mapx.get(tuple(x), x)
#    
#    class BackFrontDomain(SubDomain):
#        def inside(self, x, on_boundary):
#            # find v1x
#            return tuple(x) in problem.targetsy
#    
#    
#        # field component 1
#        def map(self, x, y):
#            y[:] = problem.vert_mapy.get(tuple(x), x)
#    
#    
#    class TopBottomDomain(SubDomain):
#        def inside(self, x, on_boundary):
#            # find v1x
#            return tuple(x) in problem.targetsz
#    
#    
#        # field component 1
#        def map(self, x, y):
#            y[:] = problem.vert_mapz.get(tuple(x), x)
#    
#    # Sub domain for Dirichlet boundary condition
#    class CenterDomain(SubDomain):
#        def inside(self, x, in_boundary):
#            return all(near(x[i], problem.center_coord[i], EPS) for i in range(problem.nDims))
  

  # Create Dirichlet boundary condition
    bcs = []
    #PKHfixed_center = DirichletBC(problem.V, Constant((0,0,0)), CenterDomain(), "pointwise")
    centerDomain = self.CenterDomain()
    centerDomain.problem = self.problem
    fixed_center = DirichletBC(problem.V, Constant((0,0,0)), centerDomain, "pointwise")
    bcs.append(fixed_center)
  
    #PKHbc1 = PeriodicBC(problem.V.sub(0), LeftRightBoundary())
    leftRightBoundary=self.LeftRightBoundary()
    leftRightBoundary.problem = self.problem
    bc1 = PeriodicBC(problem.V.sub(0), leftRightBoundary)
    bcs.append(bc1)
    #PKHbc2 = PeriodicBC(problem.V.sub(1), BackFrontDomain())
    backFrontDomain=self.BackFrontDomain()
    backFrontDomain.problem = self.problem
    bc2 = PeriodicBC(problem.V.sub(1), backFrontDomain)
    bcs.append(bc2)
    #PKHbc3 = PeriodicBC(problem.V.sub(2), TopBottomDomain())
    topBottomDomain=self.TopBottomDomain()
    topBottomDomain.problem = self.problem
    bc3 = PeriodicBC(problem.V.sub(2), topBottomDomain)
    bcs.append(bc3)
    
    problem.bcs = bcs
  
  
  
