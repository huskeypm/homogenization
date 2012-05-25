

class params:
  def __init__(self,d=1):
    self.d = d 
    self.concInitial = 0.1 # [uM]
    self.ANG_TO_UM = 1e-4
    self.dt = 0.001
    self.tStep = 20

from params import *
parms = params()

