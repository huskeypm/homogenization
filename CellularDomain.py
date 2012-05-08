# Defines boundary conditions, etc for microdomain 
from dolfin import * 

class empty:pass

class CellularDomain: 
  def __init__(self,fileMesh,fileSubdomains):
    self.problem = empty()
    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.init = 1    
    problem.name = "Cellular"

  def Setup(self):
    # mesh 
    problem = self.problem
    problem.mesh = Mesh(problem.fileMesh)
    problem.V = FunctionSpace(problem.mesh,"CG",1)
    problem.subdomains = MeshFunction(
      "uint", problem.mesh, problem.fileSubdomains)
    #self.markers = [1,4,5]

  
    

  # bcs
  def AssignBC(self):
    problem = self.problem

    u0 = Expression("1 + x[0]*x[0] + x[1]*x[1]")
    # molec boundary 
    print "Probably don't have the right BCs yet"
    bc1 = DirichletBC(problem.V,Constant(0.),problem.subdomains,1)
    # outer boundary 
    bc2 = DirichletBC(problem.V,u0,problem.subdomains,5)
    # neum
    #bc1 = DirichletBC(self.V,Constant(0.),self.subdomains,1)
    problem.bcs = [bc1,bc2]
  
 
