"
----------------------------------
Smolfin - a numerical solver of the Smoluchowski equation for interesting geometries
Copyright (C) 2012 Peter Kekenes-Huskey, Ph.D., huskeypm@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


----------------------------------------------------------------------------
"
from makeLattice import *
from homog import *
import numpy as np

class empty:pass

n = 5
particleDiams = np.arange(n)/float(n) * 2.0
print particleDiams
tau = empty()
tau.x = np.zeros(len(particleDiams))
tau.y = np.zeros(len(particleDiams))
tau.z = np.zeros(len(particleDiams))
tauo = empty()
tauo.x = np.zeros(len(particleDiams))
tauo.y = np.zeros(len(particleDiams))
tauo.z = np.zeros(len(particleDiams))
base="lattice"
parms.d


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

  r=results.molDomUnit.problem.domv  
  tauo.x[i] = r[0]
  tauo.y[i] = r[1]
  tauo.z[i] = r[2]


import matplotlib.pyplot as plt 
#plt.plot(particleDiams,tau.x,'r')
#plt.plot(particleDiams,tau.y,'g')
#plt.plot(particleDiams,tau.z,'b')
plt.plot(particleDiams,tauo.x/parms.d,'r-.')
plt.plot(particleDiams,tauo.y/parms.d,'g-.')
plt.plot(particleDiams,tauo.z/parms.d,'b-.')
plt.xlabel('Particle diams [nm]')
plt.ylabel('D$^*$/D')
plt.legend(["$D_x$","$D_y$","$D_z$"])
f = plt.gcf()
f.savefig("beadresult.png") 



