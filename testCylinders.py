"""
For making simple sphere-centric meshes with prescribed volume
fractions
"""

import sys
sys.path.append("/home/huskeypm/sources/homogenization/")
#

import numpy as np
import sys
import homoglight as hl
def Make3DGmsh(rCylinder=.25,dBox=1.,filePrefix="test"):
  """
  Creates gmsh.geo files for 2D meshes according to vol fraction 
  (only valid for densely packed cylinders) 
  foreach geo, convert to twoD mesh
               dolfin-convert volFrac_0.27.msh  volFrac_0.27.xml
  """
  lc = 0.10
  nLayers = dBox/lc  
  center = np.zeros(3)
  rBox = dBox/2. 
  cont="Point(1) = {%f,%f,%f,%f};\n" % (center[0],center[1],center[2],lc)
  cont+="Point(2) = {%f, 0, 0,%f};\n" % (rCylinder,lc)
  cont+="Point(3) = {%f, 0, 0,%f};\n" % (-rCylinder,lc)
  cont+="Point(4) = {0, %f, 0,%f};\n" % (rCylinder,lc)
  cont+="Point(5) = {0, %f, 0,%f};\n" % (-rCylinder,lc)
  cont+="Point(7) = {%f, %f, 0,%f};\n" % (rBox,rBox,lc)
  cont+="Point(8) = {%f, %f, 0,%f};\n" % (rBox,-rBox,lc)
  cont+="Point(9) = {%f, %f, 0,%f};\n" % (-rBox,rBox,lc)
  cont+="Point(10) = {%f, %f, 0,%f};\n" % (-rBox,-rBox,lc)
  cont+="""
  Circle(5) = {4, 1, 2};
  Circle(6) = {2, 1, 5};
  Circle(7) = {5, 1, 3};
  Circle(8) = {3, 1, 4};
  Line(9) = {7, 8};
  Line(10) = {8, 10};
  Line(11) = {10, 9};
  Line(12) = {9, 7};
  Line Loop(13) = {12, 9, 10, 11};
  Line Loop(14) = {5, 6, 7, 8};
  Plane Surface(15) = {13, 14};
  """
  cont+="Extrude {0,0,%f} {Surface{15}; Layers{%d};}\n" %(dBox,nLayers)
  
  name = filePrefix+".geo"
  textFile = open(name,"w")             
  textFile.write(cont)
  textFile.close()        



import os 
def computeDeffCylinder(rCylinder=0.25):
  dBox = 1.

  # create 3d  mesh 
  filePrefix="test"
  Make3DGmsh(rCylinder=rCylinder,dBox=dBox,filePrefix=filePrefix)

  # convert to msh 
  cmd = "gmsh -3 %s.geo" % filePrefix
  os.system(cmd)

  # convert to dolfin msh 
  cmd = "dolfin-convert %s.msh %s.xml" % (filePrefix,filePrefix)
  os.system(cmd) 

  # call homog light
  problem = hl.runHomog(filePrefix+".xml")
  volFrac = problem.volume/problem.volUnitCell
  deffT = problem.d_eff[0]
  deffL = problem.d_eff[2]

  return (volFrac,deffT,deffL)

def doit():
  rs = np.linspace(0.05,0.40,19)

  vfracs=np.zeros(np.shape(rs)[0])
  deffTs=np.zeros(np.shape(rs)[0])
  deffLs=np.zeros(np.shape(rs)[0])

  for i,r in enumerate(rs):
    vfrac,deffT,deffL = computeDeffCylinder(rCylinder=r)
    print deffT
    print deffL
    vfracs[i]=vfrac
    deffTs[i]=deffT
    deffLs[i]=deffL

  # plot 
  D=1.
  phi=vfracs
  hsLimit = D*phi/(2-phi)
  import matplotlib.pylab as plt
  plt.plot(phi,hsLimit,'k-',linewidth=1,label="HS limit")
  plt.plot(phi,deffTs,'k.',linewidth=2,label="Deff transverse")
  plt.plot(phi,deffLs,'k--',linewidth=1,label="Deff longitudinal")
  plt.plot(phi,phi*D,'k^',linewidth=2,label="Deff trans, analytical")
  plt.ylabel("Deff/D")
  plt.xlabel("Phi (porosity)")
  plt.title("Deff longitudinal and tranverse for cylinders") 
  plt.legend(loc=0)
  plt.gcf().savefig("cylinder.png") 
  
  

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




  doit()


