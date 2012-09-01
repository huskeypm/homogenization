# Defines boundary conditions, etc for microdomain
from dolfin import *
from params import *
from Domain import *
from util import * 


markerInsideBoundary= 1
markerOutsideBoundary= 5

class CellularUnitDomain(Domain):
  def __init__(self,fileMesh,fileSubdomains,type,\
               outpath="./",name="Cellular"):
    super(CellularUnitDomain,self).__init__(type)

    problem = self.problem
    problem.fileMesh = fileMesh
    problem.fileSubdomains = fileSubdomains
    problem.name = name 
    problem.outpath = outpath

  def Setup(self):
    # mesh
    problem = self.problem
    mesh = Mesh(problem.fileMesh)
    problem.mesh = mesh 

    # mesh is in A, so converting to um
    problem.mesh.coordinates()[:] *= parms.ANG_TO_UM

    utilObj = util(problem)
    utilObj.GeometryInitializations()
    utilObj.DefinePBCMappings()



    if(self.type=="scalar"):
        problem.V = FunctionSpace(problem.mesh,"CG",1)
    elif(self.type=="field"):
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

    print "Need to remove periodic BC"
    quit()

    bcs = []
    #PKHfixed_center = DirichletBC(problem.V, Constant((0,0,0)), CenterDomain(), "pointwise")
    centerDomain = self.CenterDomain()
    centerDomain.problem = self.problem
    fixed_center = DirichletBC(problem.V, Constant((0,0,0)), centerDomain, "pointwise")
    bcs.append(fixed_center)

    #PKHbc1 = PeriodicBC(problem.V.sub(0), PeriodicLeftRightBoundary())
    leftRightBoundary=self.PeriodicLeftRightBoundary()
    leftRightBoundary.problem = self.problem
    bc1 = PeriodicBC(problem.V.sub(0), leftRightBoundary)
    bcs.append(bc1)
    #PKHbc2 = PeriodicBC(problem.V.sub(1), PeriodicBackFrontDomain())
    backFrontDomain=self.PeriodicBackFrontDomain()
    backFrontDomain.problem = self.problem
    bc2 = PeriodicBC(problem.V.sub(1), backFrontDomain)
    bcs.append(bc2)
    #PKHbc3 = PeriodicBC(problem.V.sub(2), PeriodicTopBottomDomain())
    topBottomDomain=self.PeriodicTopBottomDomain()
    topBottomDomain.problem = self.problem
    bc3 = PeriodicBC(problem.V.sub(2), topBottomDomain)
    bcs.append(bc3)

    problem.bcs = bcs


