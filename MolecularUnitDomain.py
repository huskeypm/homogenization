# Defines boundary conditions, etc for microdomain
from dolfin import *
from params import *
from Domain import *

markerActiveSite = 1
markerMolecularBoundary =4
markerOuterBoundary=5

class MolecularUnitDomain(Domain):
  def __init__(self,fileMesh,fileSubdomains,type):
    super(MolecularUnitDomain,self).__init__(type)

    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.init = 1
    print "Enforcing gamer==1"
    self.gamer = 1
    problem.name = "Molecular"
    self.markerOuterBoundary = markerOuterBoundary

  def Setup(self):
    # mesh
    problem = self.problem
    mesh = Mesh(problem.fileMesh)
    problem.mesh = mesh
    # mesh is in A, so converting to um
    mesh.coordinates()[:] *= parms.ANG_TO_UM


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


    # geom
    self.CalcGeom(problem)


  # bcs
  def AssignBC(self,uBoundary=0):
    problem = self.problem

    print "Probably don't have the right BCs yet"
    # TODO: might need to handle BC at active site as a mixed boundary
    # condition to take into account saturation
    if(self.type=="scalar"):
        u0 = Constant(0.)
        u1 = Constant(1.)
    elif(self.type=="field"):
        u0 = Constant((0.,0,0.))
        u1 = Constant((1.,1.,1.))

    # use user-provided BC instead  
    if(uBoundary != 0):
      u1 = uBoundary


    # molec boundary
    bc0 = DirichletBC(problem.V,u0,problem.subdomains,markerActiveSite)
    # outer boundary
    bc1 = DirichletBC(problem.V,u1,problem.subdomains,markerOuterBoundary)
    problem.bcs = [bc0,bc1]


