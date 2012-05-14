# Defines boundary conditions, etc for cell macrodomain
from dolfin import *
from params import *
from Domain import *

class CellularDomain(Domain):
  def __init__(self,fileMesh,fileSubdomains):
    super(CellularDomain,self).__init__()

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
    if(self.type=="scalar"):
        problem.V = FunctionSpace(problem.mesh,"CG",1)
    elif(self.type=="field"):
        problem.V = VectorFunctionSpace(problem.mesh,"CG",1)

    problem.subdomains = MeshFunction(
      "uint", problem.mesh, problem.fileSubdomains)
    #self.markers = [1,4,5]




  # bcs
  def AssignBC(self):
    problem = self.problem

    print "Add in LCC"
    #if(self.type=="scalar"):
    #    dudn = Constant(2.)
    #elif(self.type=="field"):
    #    dudn = Constant((2.,2.,2.))
    dudn = Constant(2.)

    # Currently don't have any Dirichlet BCs
    problem.bcs=[]

    problem.dudn = dudn  # assigning a flux on the entire boundary for now

