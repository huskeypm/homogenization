# Defines boundary conditions, etc for microdomain
from dolfin import *
from params import *
from Domain import *

class CellularUnitDomain(Domain):
  def __init__(self,fileMesh,fileSubdomains):
    super(CellularUnitDomain,self).__init__()

    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.init = 1
    problem.name = "Cellular"

  def Setup(self,type="scalar"):
    # mesh
    problem = self.problem
    problem.mesh = Mesh(problem.fileMesh)
    self.type = type
    if(type=="scalar"):
        problem.V = FunctionSpace(problem.mesh,"CG",1)
    elif(type=="field"):
        problem.V = VectorFunctionSpace(problem.mesh,"CG",1)

    problem.subdomains = MeshFunction(
      "uint", problem.mesh, problem.fileSubdomains)
    #self.markers = [1,4,5]

    # geom
    self.CalcGeom(problem)




  # bcs
  def AssignBC(self):
    problem = self.problem

    if(self.type=="scalar"):
        u0 = Constant(0.)
        #u1 = Expression("1 + x[0]*x[0] + x[1]*x[1]")
        # I don't think I want any specific BCs here for homog prob.
        u1 = Constant(1.)
    elif(self.type=="field"):
        u0 = Constant((0.,0,0.))
        u1 = Constant((1.,1.,1.))


    print "Not sure if I should be using any Dirichlet BC"
    bc0 = DirichletBC(problem.V,u0,problem.subdomains,parms.markerActiveSite)
    bc1 = DirichletBC(problem.V,u1,problem.subdomains,parms.markerOuterBoundary)
    # neum

    problem.bcs = [bc0,bc1]


