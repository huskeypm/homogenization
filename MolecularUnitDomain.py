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

#class TestRL(SubDomain):
#  def inside(self, x, on_boundary):
#    cond = (np.abs(x[0]- -1) < EPS or np.abs(x[0]-1) < EPS)
#    #if(on_boundary):
#    #  print "x",x[0],cond
#    return (on_boundary and cond)
#
#class TestTB(SubDomain):
#  def inside(self, x, on_boundary):
#    cond = (np.abs(x[1]- -1) < EPS or np.abs(x[1]-1) < EPS)
#    #if(on_boundary):
#    #  print "y",x[1],cond
#    return (on_boundary and cond)
#
#class TestBF(SubDomain):
#  def inside(self, x, on_boundary):
#    cond = (np.abs(x[2]- -1) < EPS or np.abs(x[2]-1) < EPS)
#    #if(on_boundary):
#    #  print "z",x[2],cond
#    return (on_boundary and cond)


#
# filePotential - electrostatic potential from APBS, interpolated to FE mesh
class MolecularUnitDomain(Domain):
  def __init__(self,fileMesh,fileSubdomains,filePotential="none",\
    # scalar (wrong) or field
    type="field",\
    # doe mesh come from gamer?
    gamer=1,\
    outpath="./",\
    name = "Molecular",\
    ):
    super(MolecularUnitDomain,self).__init__(type)

    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.filePotential= filePotential 
    problem.init = 1
    #print "Enforcing gamer==1"
    self.gamer = gamer
    problem.name = name
    problem.outpath = outpath
    self.markerOuterBoundary = markerOuterBoundary

    if(problem.filePotential!="none"):
      problem.smolMode="true"
    else:
      problem.smolMode="false"

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
    if(problem.smolMode=="true"):       
      print "PUT IN SEPARATE FUNCTION"
      Vtemp = FunctionSpace(mesh,"CG",1)
      psi = Function(Vtemp,problem.filePotential)
      smol.ElectrostaticPMF(problem,psi,V=Vtemp,q=q) # pmf stored internally             
    
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
    centerDomain = self.CenterDomain()
    centerDomain.problem = self.problem
    fixed_center = DirichletBC(problem.V, Constant((0,0,0)), centerDomain, "pointwise")
    #bcs.append(fixed_center)

    #PKH 120901 leftRightBoundary=self.PeriodicLeftRightBoundary()
    leftRightBoundary=self.LeftRightBoundary()
    leftRightBoundary.problem = self.problem
    #PKH 120901 bc1 = PeriodicBC(problem.V.sub(0), leftRightBoundary)
    tc1 = DirichletBC(problem.V.sub(0), Constant(1.),leftRightBoundary)
    bc1 = DirichletBC(problem.V.sub(0), Constant(0.),leftRightBoundary)
    bcs.append(bc1)

    #PKH 120901 backFrontBoundary=self.PeriodicBackFrontBoundary()
    backFrontBoundary=self.BackFrontBoundary()
    backFrontBoundary.problem = self.problem
    #PKH 120901 bc2 = PeriodicBC(problem.V.sub(1), backFrontBoundary)
    tc2 = DirichletBC(problem.V.sub(1), Constant(1.),backFrontBoundary)
    bc2 = DirichletBC(problem.V.sub(1), Constant(0.),backFrontBoundary)
    bcs.append(bc2)

    #PKH 120901 topBottomBoundary=self.PeriodicTopBottomBoundary()
    topBottomBoundary=self.TopBottomBoundary()
    topBottomBoundary.problem = self.problem
    #PKH 120901 bc3 = PeriodicBC(problem.V.sub(2), topBottomBoundary)
    tc3 = DirichletBC(problem.V.sub(2), Constant(1.),topBottomBoundary)
    bc3 = DirichletBC(problem.V.sub(2), Constant(0.),topBottomBoundary)
    bcs.append(bc3)

    testBC=1
    if(testBC==1):
      z = Function(problem.V)
      tc1.apply(z.vector())
      tc2.apply(z.vector())
      tc3.apply(z.vector())
      File("appliedBCs.pvd") << z
    
    problem.bcs = bcs
  
  
  
