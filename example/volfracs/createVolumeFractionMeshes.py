"""
For making simple sphere-centric meshes with prescribed volume
fractions
"""

import sys
#sys.path.append("/home/huskeypm/sources/smolfin/")
from gamer import SurfaceMesh, GemMesh
from meshmanipulate import *
from gmshUtil import *

#
# Revisions
#       10.08.10 inception
#

import numpy as np
def Make3DMesh(volFrac=0.27,rSphere=12.5,dBox=65.,name="test"):
  sm = SurfaceMesh("./meshes/cube_high.off")
  center(sm)

  rBox = dBox/2.
  scaleGeom(sm,rBox)

  inner = SurfaceMesh(4)
  scaleGeom(inner,rSphere)

  inner.as_hole = True
  gem_mesh = GemMesh([sm, inner])
  gem_mesh.write_dolfin(name)

import sys



def doit(fileIn):
  rSphere = 12.5

  volFracs = [0.10,.20,.27,.34,.50]

  for i,volFrac in enumerate(volFracs):
    vSphere = (4/3.*np.pi * rSphere**3)
    vBox = vSphere/volFrac
    dBox = vBox**(1/3.)
    #name = "volFrac_%4.2f.geo" % volFrac
    # test with 2D MakeGmsh(volFrac=volFrac,rSphere=rSphere,dBox=dBox)
    name = "./meshes/volFrac_%4.2f.xml" % volFrac
    Make3DMesh(volFrac=volFrac,rSphere=rSphere,dBox=dBox,name=name)
    



if __name__ == "__main__":
  msg="""
Purpose: 
  Creates volumetric meshes for different vol fractions  

Usage:
  script.py run   

Notes:

"""
  remap = "none"


  import sys
  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-arg1"):
      arg1=sys.argv[i+1] 




  doit(fileIn)


