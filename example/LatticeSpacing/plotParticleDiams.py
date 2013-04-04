

import numpy as np
import matplotlib.pyplot as plt

def doit(modes):

  for i,mode in enumerate(modes): 
    d = np.loadtxt(mode+"/outs")
    if(mode=="skeletal"):
      plt.plot(d[:,0],d[:,1],"k-",label="Skeletal ($d_{MA}=16.5 \mu m$)")
    if(mode=="cardiac"):
      plt.plot(d[:,0],d[:,1],"k--",label="Cardiac ($d_{MA}=18.3 \mu m$)")
  
  plt.title("Tortuosity")
  plt.ylabel("$D*/D_{Bulk}$")
  plt.xlabel("Particle diameter [$nm$]")
  plt.legend()
  plt.ylim([0,1.0])
  plt.xlim([0,9.5])
  
  c = plt.gcf()
  c.savefig("unitcell_particlediams.png")


import sys
#
# Revisions
#       10.08.10 inception
#

if __name__ == "__main__":
  msg="""
Purpose: 

Usage:
  script.py <arg>

Notes:
"""
  remap = "none"


  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  modes = []
  for i,arg in enumerate(sys.argv):
    if(arg=="-skeletal"):
      modes.append("skeletal")
    if(arg=="-cardiac"):
      modes.append("cardiac")




  doit(modes)


