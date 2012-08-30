# Defines boundary conditions, etc for cell macrodomain
from dolfin import *
from params import *
from Domain import *

class CellularDomain(Domain):
  def __init__(self,fileMesh,fileSubdomains,type):
    super(CellularDomain,self).__init__(type)

    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.name = "Cellular"

  def Setup(self):
    # mesh
    problem = self.problem
    problem.mesh = Mesh(problem.fileMesh)

    if(self.type=="scalar"):
        problem.V = FunctionSpace(problem.mesh,"CG",1)
    elif(self.type=="field"):
        problem.V = VectorFunctionSpace(problem.mesh,"CG",1)

    problem.subdomains = MeshFunction(
      "uint", problem.mesh, problem.fileSubdomains)


    # geom
    self.CalcGeom(problem)

  # bcs
  def AssignBC(self):
    problem = self.problem

    print "Add in LCC"
    dudn = Constant(2.)

    # Currently don't have any Dirichlet BCs
    problem.bcs=[]

    #problem.dudn = dudn  # assigning a flux on the entire boundary for now
    t = 0
    problem.dudn = Expression("-2 * exp(-t)",
      t = t)
    #problem.dudn = Constant(2.)

