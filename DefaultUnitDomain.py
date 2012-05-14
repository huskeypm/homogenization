# Defines boundary conditions, etc for microdomain
from dolfin import *
from params import *
from Domain import *

# bcs
def boundary(x,on_boundary):
  return on_boundary

class DefaultUnitDomain(Domain):
  def __init__(self):
    super(DefaultUnitDomain,self).__init__()

    problem = self.problem
    problem.init = 1
    problem.name = "Default"

  def Setup(self,type="scalar"):
    # mesh
    problem = self.problem
    problem.mesh = UnitCube(8,8,8)
    self.type = type
    if(self.type=="scalar"):
        problem.V = FunctionSpace(problem.mesh,"CG",1)
    elif(self.type=="field"):
        problem.V = VectorFunctionSpace(problem.mesh,"CG",1)

    # geom
    self.CalcGeom(problem)



  def AssignBC(self):
    problem = self.problem

    print "Probably don't have the right BCs yet"
    if(self.type=="scalar"):
        u1 = Constant(1.)
    elif(self.type=="field"):
        u1 = Constant((1.,1.,1.))

    bc = DirichletBC(problem.V, u1, boundary)
    problem.bcs = [bc]


