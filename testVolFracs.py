import sys
import numpy as np
import matplotlib.pylab as plt
#sys.path.append("/home/huskeypm/sources/smolfin/")
#sys.path.append("/home/huskeypm/sources/homogenization/")

from dolfin import *
#import imp
#pb  = imp.load_source('poissonboltzmann.py', '/home/huskeypm/sources/smolfin/')
import poissonboltzmann as pb
meshName = "meshes/volFrac_0.10.xml"

# quick test
#mesh = Mesh(meshName)
#import poissonboltzmann as pb
#pb.params.molRad = 12.5
#pb.params.domRad = 200
#(V,x) = pb.SolvePoissonBoltzmann(mesh,meshType="gamer")

import homoglight as hl
debug = 0
zProt = -3 # [mV]
molRad = 12.5 # [A]5 # [A] 

#import homog_w_electro as hwe
  #Dx = hwe.doit(meshFile="meshes/volFrac_0.10.xml",meshType="gamer")

# molRad - radius of enzeyme 
# z - enzyme charge
# q - substrate charge 

def CalcPMF(mesh,molRad,zLig,zProt,meshType="dolfin"):
   print "zLig=",zLig
   print "zProt",zProt
   pb.parms.zLig = zLig
   pb.parms.zProt = zProt
   pb.parms.molRad = 12.5
   pb.parms.sphericalDomainBoundary=False
   pb.parms.update()
   (Vdumm,psi) = pb.SolvePoissonBoltzmann(mesh,meshType=meshType)
   name = "elect.pvd"
   psiar = np.asarray(psi.vector()[:])
   print "Psi min/max %f/%f [mV?]" %(np.min(psiar),np.max(psiar))
    
   File(name) << psi

   # convert from [mV] to [kT/|z|]  
   pmf = Function(Vdumm)
   pmf.vector()[:]= psi.vector()[:]*pb.parms.F_o_RT*pb.parms.kT
   pmfar = np.asarray(pmf.vector()[:])
   print "Pmf min/max %f/%f [mV?]" %(np.min(pmfar),np.max(pmfar))
   return pmf 
  

def call(meshPrefix,zLig,volFrac,molGamer=0):
#  hwe.params.z = -1
#  hwe.pb.molRad = 12.5
#  hwe.pb.sphericalDomainBoundary=False
#  hwe.params.molRad = 12.5
#  hwe.params.sphericalDomainBoundary=False
#  hwe.params.q=1.


  # mesh 
  meshName = meshPrefix+"_mesh.xml.gz"
  mesh = Mesh(meshName)

  # ESP 
  #psi = CalcPMF(mesh,molRad,z,q,meshType="gamer")
  psi = CalcPMF(mesh,molRad,zLig,zProt)

  
  # double check size just in case
  mn = np.min(mesh.coordinates(),axis=0)
  mx = np.max(mesh.coordinates(),axis=0)
  boxVol = np.prod(mx-mn)
  sphereVol = 4/3.*np.pi*(molRad**3)
  print meshName
  print "Vol frac %f vs %f " % (sphereVol/boxVol,volFrac)


  #molDomUnit.smolMode = smolMode
  #solve_homogeneous_unit(molDomUnit,type="field",debug=debug,smolMode=smolMode)
  #homog.SolveHomogSystem(molPrefix=meshPrefix,molGamer=molGamer)  
  #quit()
  #results = homog.SolveHomogSystem(molPrefix=meshPrefix,molGamer=molGamer,\
  #  smolMode=True,smolPsi=pmf,smolq = zLig)
  results = hl.runHomog(fileXML=meshName,psi=psi,q=zLig,smolMode=True)
  Dx = results.d_eff[0]

  #(V,x) = pb.SolvePoissonBoltzmann(mesh)
  # do electro homog
  # store value
  #results[i,j] = np.random.rand() +  i*j

  # moleculardomain.q =q 
  # moleculardomain.psi = psi
  
  #Dx = hwe.doit(meshFile=meshName, meshType="gamer")
  return Dx
  
  
def validation():
  meshPrefix = "example/volfracs/volFrac_0.10"       
  volFrac = 0.1
  zLig=1
  zLig=-1
  valuep =  call(meshPrefix,zLig,volFrac,molGamer=0)             
  quit()
  zLig=0
  value0 =  call(meshPrefix,zLig,volFrac,molGamer=0)             
  zLig=-1
  valuen =  call(meshPrefix,zLig,volFrac,molGamer=0)             
  #value130404=0.535335
  value130520=0.55234253        

  print "sphere, neutral ", value0
  print "sphere, positive ", valuep
  print "sphere, negative ", valuen

  assert(np.abs(valuep-value130520) < 0.001), "RESULT CHANGED. DO NOT COMMIT"




def doit(asdf):
  ## params 

  volFracs = np.array([0.1,0.2,0.27,0.34,0.5])
  #meshes = np.array([0.1,0.2,0.27,0.34,0.5])
  meshes = volFracs
  # kappa hard coded into PB solver
  qs = np.array([-2,-1,0,1,2])
  #qs = np.array([-1,0])           
  
  
  if(debug==1):
    qs=[0]
    volFracs=[0.1]
    #qs=[2]
    volFracs=[0.27]
    meshes = volFracs

  results = np.zeros([ np.shape(volFracs)[0], np.shape(qs)[0] ])
  
  for i, volFrac in enumerate(volFracs):
    for j, q in enumerate(qs):
      #meshName = "meshes/volFrac_%4.2f.xml" % meshes[i]
      #meshName = "meshes/volFrac_%4.2f_mesh.xml.gz" % meshes[i]
      meshPrefix= "example/volfracs/volFrac_%4.2f" % meshes[i]
      #molGamer = 1
      molGamer = 0
      results[i,j] = call(meshPrefix,q,volFrac,molGamer=molGamer)

  print "WARNING: should obtain these directly from sims, not ests"
  volFracActual = 1-volFracs

  #if(debug==1):
  #  return
  
  
  plt.figure()
  col = ["r--","r-","k-","b-","b--"]
  for j, q in enumerate(qs):
    label = "zLig= %d " % q
    plt.plot(volFracActual,results[:,j], col[j],label=label)
  
  title="Effective Diffusion versus volume fraction (zProtein=%5.2f)"\
     % zProt
  plt.title(title)
  plt.xlabel("Volume fraction")
  plt.ylabel("Dx/Dbulk")
  plt.legend(loc=0)
  
  plt.gcf().savefig("volfrac.png")
  








import sys
#
# Revisions
#       10.08.10 inception
#


if __name__ == "__main__":
  import sys
  scriptName= sys.argv[0]
  msg="""
Purpose: 
  Compute diffusion constant for a variety of geometries 
 
Usage:
"""
  msg+="  %s -debug/-validation" % (scriptName)
  msg+="""
  
 
Notes:

"""
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  for i,arg in enumerate(sys.argv):
    if(arg=="-debug"):
      debug=1
    if(arg=="-validation"):
      validation()
      quit()




  doit(fileIn)


