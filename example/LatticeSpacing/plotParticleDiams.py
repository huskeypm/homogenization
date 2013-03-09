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
import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("outs")
plt.plot(d[:,0],d[:,1],"k-",label="Rest (distMA=16.5 um)")

plt.title("Tortuosity")
plt.xlabel("D* [$um^2$/s]")
plt.ylabel("D* [$um^2$/s]")
plt.xlabel("Particle diameter [$nm$]")
plt.legend()
plt.ylim([0,1.0])
plt.xlim([0,8.0])

c = plt.gcf()
c.savefig("unitcell_particlediams.png")


