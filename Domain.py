from dolfin import *
class empty:pass

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



