# Defines boundary conditions, etc for microdomain
from dolfin import *

class empty:pass

class MolecularUnitDomain:
  def __init__(self,fileMesh,fileSubdomains):
    self.problem = empty()
    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.init = 1
    problem.name = "Molecular"

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

    print "Probably don't have the right BCs yet"
    if(self.type=="scalar"):
        u0 = Constant(0.)
        u1 = Constant(1.)
    elif(self.type=="field"):
        u0 = Constant((0.,0,0.))
        u1 = Constant((1.,1.,1.))

    # molec boundary
    bc0 = DirichletBC(problem.V,u0,problem.subdomains,1)
    # outer boundary
    bc1 = DirichletBC(problem.V,u1,problem.subdomains,5)
    # neum
    #bc1 = DirichletBC(self.V,Constant(0.),self.subdomains,1)
    problem.bcs = [bc0,bc1]


