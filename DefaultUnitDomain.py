# Defines boundary conditions, etc for microdomain
from dolfin import *
from params import *
from Domain import *

# bcs
def boundary(x,on_boundary):
  return on_boundary

class DefaultUnitDomain(Domain):
  def __init__(self,type):
    super(DefaultUnitDomain,self).__init__(type)

    problem = self.problem
    problem.name = "Default"

  def Setup(self):
    # mesh
    problem = self.problem
    problem.mesh = UnitCube(8,8,8)
    if(self.type=="scalar"):
        problem.V = FunctionSpace(problem.mesh,"CG",1)
    elif(self.type=="field"):
        problem.V = VectorFunctionSpace(problem.mesh,"CG",1)

    # geom
    self.CalcGeom(problem)



  # Have option of passing in uBoundary 
  def AssignBC(self,uBoundary=0):
    problem = self.problem

    #print "Probably don't have the right BCs yet"
    if(self.type=="scalar"):
        u1 = Constant(1.)
    elif(self.type=="field"):
        u1 = Constant((1.,1.,1.))

    # use user-provided BC instead  
    if(uBoundary != 0):
      u1 = uBoundary

    bc = DirichletBC(problem.V, u1, boundary)
    problem.bcs = [bc]

    # Flux condition (just a placeholder for debugging)
    t = 0
    problem.dudn = Expression("0.5  * exp(-t/2.)",
      t=t)


