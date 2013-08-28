# very simple module for homogenzation 

from homog import *
def runHomog(fileXML="test.xml",psi="none",smolMode=False,q=0,verbose=False,gamer=0):
  fileSubdomains = "none"
  molDomUnit = MolecularUnitDomain(fileXML,fileSubdomains,gamer=gamer,\
    psi=psi,q=q)
  molDomUnit.Setup()
  molDomUnit.AssignBC()

  solve_homogeneous_unit(molDomUnit,smolMode=smolMode)

  problem = molDomUnit.problem
  if(verbose):
    print problem.volume
    print problem.volUnitCell
    print problem.d_eff
  return problem 



#!/usr/bin/env python
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
  Light-weight wrapper for running homogeniation on xml files 
 
Usage:
"""
  msg+="  %s -file <filename>" % (scriptName)
  msg+="""
  
 
Notes:

"""
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  fileIn= sys.argv[1]
  if(len(sys.argv)==3):
    print "arg"

  fileXML = "none" 
  for i,arg in enumerate(sys.argv):
    if(arg=="-file"):
      fileXML=sys.argv[i+1] 



  runHomog(fileXML,verbose=True)



