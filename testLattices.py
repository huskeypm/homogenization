from makeLattice import *
from homog import *
import numpy as np

class empty:pass

n = 1
particleDiams = np.arange(n)/float(n) * 2.0
print particleDiams
tau = empty()
tau.x = np.zeros(len(particleDiams))
tau.y = np.zeros(len(particleDiams))
tau.z = np.zeros(len(particleDiams))
base="lattice"


i=0
for i in range(n):         
  name = base+"%.2d" % (i)
  print "### Running" + name 
  MakeDomain(base=name,scale=1.0,particleDiam=particleDiams[i])
  results = SolveHomogSystem(debug=0,\
      root="./",\
      cellPrefix="none", molPrefix=name,wholeCellPrefix="none",\
      #smolMode = False,
      molGamer=1)
  r=results.molDomUnit.problem.d_eff
  r = r/np.linalg.norm(r)
  tau.x[i] = r[0]
  tau.y[i] = r[1]
  tau.z[i] = r[2]

import matplotlib.pyplot as plt 
plt.plot(particleDiams,tau.x,'r')
plt.plot(particleDiams,tau.y,'g')
plt.plot(particleDiams,tau.z,'b')
f = plt.gcf()
f.savefig("result.png") 



