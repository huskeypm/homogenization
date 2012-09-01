from dolfin import *
import numpy as np
class empty:pass

EPS = 1.e-10
# make pretty lax
EPS = 1.e-1  

class Domain(object):
  #
  # CLASSES
  #
  class LeftRightBoundary(SubDomain):
          def inside(self, x, on_boundary):
              problem = self.problem 
              return ( np.abs(x[0]-problem.boundsMin[0]) < EPS or np.abs(x[0]-problem.boundsMax[0]) < EPS)

      
  class BackFrontDomain(SubDomain):
          def inside(self, x, on_boundary):
              problem = self.problem 
              return ( np.abs(x[1]-problem.boundsMin[1]) < EPS or np.abs(x[1]-problem.boundsMax[1]) < EPS)
      
  class TopBottomDomain(SubDomain):
          def inside(self, x, on_boundary):
              problem = self.problem 
              return ( np.abs(x[2]-problem.boundsMin[2]) < EPS or np.abs(x[2]-problem.boundsMax[2]) < EPS)
      
  class PeriodicLeftRightBoundary(SubDomain):
          def inside(self, x, on_boundary):
              # find v1x
              problem = self.problem 
              return tuple(x) in problem.targetsx
      
      
          # field component 1
          def map(self, x, y):
              problem = self.problem 
              y[:] = problem.vert_mapx.get(tuple(x), x)
      
  class PeriodicBackFrontDomain(SubDomain):
          def inside(self, x, on_boundary):
              # find v1x
              problem = self.problem 
              return tuple(x) in problem.targetsy
      
      
          # field component 1
          def map(self, x, y):
              problem = self.problem 
              y[:] = problem.vert_mapy.get(tuple(x), x)
      
      
  class PeriodicTopBottomDomain(SubDomain):
          def inside(self, x, on_boundary):
              # find v1x
              problem = self.problem 
              return tuple(x) in problem.targetsz
      
      
          # field component 1
          def map(self, x, y):
              problem = self.problem 
              y[:] = problem.vert_mapz.get(tuple(x), x)
      
      # Sub domain for Dirichlet boundary condition
  class CenterDomain(SubDomain):
          def inside(self, x, in_boundary):
              problem = self.problem 
              return all(near(x[i], problem.center_coord[i], EPS) for i in range(problem.nDims))



  #
  # FUNCTIONS
  # 
  def __init__(self,type):
    problem = empty()
    problem.gamma = 1.
    problem.volUnitCell= 1.
    self.type = type
    problem.init = 1
    self.gamer = 0 # meshes recorded by gamer are marked differently

    self.problem = problem

  def AssignBC(self,uBoundary=0):
    problem = self.problem

    #print "Probably don't have the right BCs yet"
    # TODO: might need to handle BC at active site as a mixed boundary
    # condition to take into account saturation
    if(self.type=="scalar"):
        print "REPLACE THIS AS DEBUG OPTION. scalar approach is nonsensical"
        quit()
        u0 = Constant(0.)
        u1 = Constant(1.)
    elif(self.type=="field"):
        u0 = Constant((0.,0,0.))
        u1 = Constant((1.,1.,1.))

    # use user-provided BC instead  
    if(uBoundary != 0):
      u1 = uBoundary

  # Create Dirichlet boundary condition
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



  def CalcGeom(self,problem):
    # SA
    if(self.gamer==0):
      areaExpr = Constant(1.) * ds
    if(self.gamer==1):
      print "WARNING: need to automatically integtrate over all surfs"
      areaExpr = Constant(1.)*ds(1) + Constant(1.)*ds(5)
    area = assemble(areaExpr, mesh=problem.mesh)

    problem.surfaceArea = area
    # this is the 'beta' term in the Goel papers 
    print "SA: %e [um^2]" % area

    # VOL 
    if(self.gamer==0):
      vol = assemble(Constant(1.) * dx, mesh=problem.mesh)
    if(self.gamer==1):
      vol = assemble(Constant(1.) * dx(1), mesh=problem.mesh)
    problem.volume = vol
    print "Volume: %e [um^3]" % vol  

