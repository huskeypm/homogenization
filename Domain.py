from dolfin import *
class empty:pass

EPS = 1.e-10

class Domain(object):
  def __init__(self,type):
    problem = empty()
    problem.gamma = 1
    self.type = type
    problem.init = 1
    self.gamer = 0 # meshes recorded by gamer are marked differently

    self.problem = problem


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

  class LeftRightBoundary(SubDomain):
          def inside(self, x, on_boundary):
              # find v1x
              problem = self.problem 
              return tuple(x) in problem.targetsx
      
      
          # field component 1
          def map(self, x, y):
              problem = self.problem 
              y[:] = problem.vert_mapx.get(tuple(x), x)
      
  class BackFrontDomain(SubDomain):
          def inside(self, x, on_boundary):
              # find v1x
              problem = self.problem 
              return tuple(x) in problem.targetsy
      
      
          # field component 1
          def map(self, x, y):
              problem = self.problem 
              y[:] = problem.vert_mapy.get(tuple(x), x)
      
      
  class TopBottomDomain(SubDomain):
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



