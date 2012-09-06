from makeLattice import *
from homog import *
import numpy as np

particleDiams = np.arange(10)/5.
tau = np.zeros(len(particleDiams))
name="x"


i=0
for i in particleDiams:
  MakeDomain(base=name,scale=1.0,particleDiam=particleDiams[i])
  results = SolveHomogSystem(debug=0,\
      root="./",\
      cellPrefix="none", molPrefix=name,wholeCellPrefix="none",\
      #smolMode = False,
      molGamer=1)
  r=results.molDomUnit.problem.d_eff
  r = r/np.linalg.norm(r)
  tau[i] = r[0]

import matplotlib.pyplot as plt 

plt.plot(particleDiams,tau)
f = plt.gcf()
f.savefig("result.png") 



