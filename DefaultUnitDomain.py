# Defines boundary conditions, etc for microdomain
from dolfin import *
from params import *
from Domain import *
from util import *


# bcs
def boundary(x,on_boundary):
  return on_boundary

class DefaultUnitDomain(Domain):
  def __init__(self,type="field"):
    super(DefaultUnitDomain,self).__init__(type)

    problem = self.problem
    problem.name = "Default"

  def Setup(self):
    # mesh
    problem = self.problem
    #problem.mesh = UnitCube(8,8,8)
    problem.mesh = UnitSphere(8)
    if(self.type=="scalar"):
        problem.V = FunctionSpace(problem.mesh,"CG",1)
    elif(self.type=="field"):
        problem.V = VectorFunctionSpace(problem.mesh,"CG",1)

    utilObj = util(problem)
    utilObj.GeometryInitializations()
    utilObj.DefinePBCMappings()

    # geom
    self.CalcGeom(problem)


