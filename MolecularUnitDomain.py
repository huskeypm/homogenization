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


#
# filePotential - electrostatic potential from APBS, interpolated to FE mesh
class MolecularUnitDomain(Domain):
  def __init__(self,fileMesh,fileSubdomains,filePotential="none",\
    # scalar (wrong) or field
    type="field",\
    # doe mesh come from gamer?
    gamer=1\
    ):
    super(MolecularUnitDomain,self).__init__(type)

    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.filePotential= filePotential 
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
    utilObj.GeometryInitializations()
    utilObj.DefinePBCMappings()

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
    if(problem.filePotential!="none"):
      print "PUT IN SEPARATE FUNCTION"
      Vtemp = FunctionSpace(mesh,"CG",1)
      psi = Function(Vtemp,problem.filePotential)
      smol.ElectrostaticPMF(problem,Vtemp,psi,q=q) # pmf stored internally             
    
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
  
  
  
