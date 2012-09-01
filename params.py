

class params:
  def __init__(self,d=2.): # Be sure to use float 
    self.d = d 
    self.concInitial = 0.1 # [uM]
    #self.ANG_TO_UM = 1e-4
    print "WARNING: removed ANG TO UM" 
    self.ANG_TO_UM = 1     
    self.dt = 0.001
    self.tStep = 20

from params import *
parms = params()

