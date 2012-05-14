from dolfin import *
class empty:pass

class Domain(object):
  def __init__(self):
    problem = empty()
    problem.gamma = 1

    self.problem = problem


  def CalcGeom(self,problem):
    # SA
    area = assemble(Constant(1.) * ds, mesh=problem.mesh)
    problem.surfaceArea = area
    # this is the 'beta' term in the Goel papers 
    print "SA: %f" % area

    # VOL 
    vol = assemble(Constant(1.) * dx, mesh=problem.mesh)
    problem.volume = vol
    print "Volume: %f" % vol  



