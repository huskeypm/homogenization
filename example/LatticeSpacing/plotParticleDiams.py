

import numpy as np
import matplotlib.pyplot as plt

def doit():

  path = '/net/data/huskeypm/GiKe12a/130118/'
  cardiac = np.loadtxt(path+'cardiac/cardiac.out',usecols=[0,1,3])
  skeletal  = np.loadtxt(path+'skeletal/skeletal.out',usecols=[0,1,3])
  
  def process(case):
    list = []
    l=np.shape(np.where(case[:,1] == 0.1))[1]
    cases = case[np.where(case[:,1] == 0.1),0]
    list.append(case[np.where(case[:,1] == 0.1),2])
    list.append(case[np.where(case[:,1] == 0.15),2])
    list.append(case[np.where(case[:,1] == 0.2),2])
    list = np.reshape(np.asarray(list).T,[l,3])
    return cases,list 
  
  scases,skeletalList = process(skeletal)
  ccases,cardiacList = process(cardiac)
  
  r=np.ndarray.flatten(scases)
  m=np.ndarray.flatten(np.mean(skeletalList,axis=1))
  yerr=1.96*np.std(skeletalList,axis=1)/np.sqrt(3)
  plt.errorbar(r,m,yerr=yerr,fmt='ko',capsize=5,label="Skeletal ($d_{MA}=16.5 nm$)")
  
  r=np.ndarray.flatten(ccases)
  m=np.ndarray.flatten(np.mean(cardiacList,axis=1))
  yerr=1.96*np.std(cardiacList,axis=1)/np.sqrt(3)
  plt.errorbar(r,m,yerr=yerr,fmt='ro',capsize=5,label="Cardiac ($d_{MA}=18.3 nm$)")
  
  plt.title("Tortuosity")
  plt.ylabel("$D*/D_{Bulk}$")
  plt.xlabel("Particle diameter [$nm$]")
  plt.legend()
  plt.ylim([0,1.0])
  plt.xlim([0,9.5])
  
  print "printing %s " % (  path+"unitcell_particlediams.png")
  plt.gcf().savefig(path+"unitcell_particlediams.png") 
  
  
  



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


  doit()


