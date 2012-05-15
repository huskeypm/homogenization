# Defines boundary conditions, etc for microdomain
# See notes from 120514_grids.tex
from dolfin import *
from params import *
from Domain import *

def boundary(x,on_boundary):
  return on_boundary


class CellularUnitDomain_TnC(Domain):
  def __init__(self):
    super(CellularUnitDomain_TnC,self).__init__()

    problem = self.problem
    problem.init = 1
    problem.name = "Cellular"

  def Setup(self,type="scalar"):
    # mesh
    problem = self.problem
    problem.mesh = UnitCube(8,8,8)

    # resizing mesh to be 0.05x0.05x0.05 [um]
    problem.mesh.coordinates()[:] *= 0.05/2. # since centered at origin


    self.type = type
    if(type=="scalar"):
        problem.V = FunctionSpace(problem.mesh,"CG",1)
    elif(type=="field"):
        problem.V = VectorFunctionSpace(problem.mesh,"CG",1)

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
    bc1 = DirichletBC(problem.V,u1,boundary)
    # neum

    problem.bcs = [bc1]


