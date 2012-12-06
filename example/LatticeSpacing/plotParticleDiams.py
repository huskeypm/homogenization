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


