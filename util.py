from dolfin import *
import numpy as np

###
### Class for dealing with various geometric issues (bounds, pairing between vertices) 
### 

class util: 
  def __init__(self,problem):
    self.prob = problem 

  def CalcBounds(self,mesh):
    prob = self.prob
  
    coords = mesh.coordinates()
    prob.nDims = (np.shape(coords))[1]
  
    boundMin=np.zeros(prob.nDims)
    boundMax=np.zeros(prob.nDims)
    boundIdxMin=np.zeros(prob.nDims,"int")
    boundIdxMax=np.zeros(prob.nDims,"int")
  
    boundMin[0] = coords[:,0].min()
    boundIdxMin[0] = coords[:,0].argmin()
    boundMax[0] = coords[:,0].max()
    boundIdxMax[0] = coords[:,0].argmax()
    
    boundMin[1] = coords[:,1].min()
    boundIdxMin[1] = coords[:,1].argmin()
    boundMax[1] = coords[:,1].max()
    boundIdxMax[1] = coords[:,1].argmax()
    
    boundMin[2] = coords[:,2].min()
    boundIdxMin[2] = coords[:,2].argmin()
    boundMax[2] = coords[:,2].max()
    boundIdxMax[2] = coords[:,2].argmax()
  
    return (boundMin,boundIdxMin,boundMax,boundIdxMax)
  
  def CalcMidpoint(self,mesh):
    (boundMin,boundIdxMin,boundMax,boundIdxMax) = self.CalcBounds(mesh)
    return (boundMin + boundMax)/2.
  
  def CalcRanges(self,mesh):
    (boundMin,boundIdxMin,boundMax,boundIdxMax) = self.CalcBounds(mesh)
    prob = self.prob
  
   
    ranges = np.zeros(prob.nDims)
    ranges[0] = boundMax[0]-boundMin[0]
    ranges[1] = boundMax[1]-boundMin[1]
    ranges[2] = boundMax[2]-boundMin[2]
  
    return ranges
  
  def CenterMesh(self,mesh):
  
    (boundMin,boundIdxMin,boundMax,boundIdxMax) = self.CalcBounds(mesh)
    mp = self.CalcMidpoint(mesh)
    #print CalcMidpoint(mesh)

    center = 1
    if(center):
      for i,c in enumerate( mesh.coordinates() ):
        c -= mp
        mesh.coordinates()[i] = c
  
  def DefinePBCMappings(self):
    prob = self.prob 
    mesh = prob.mesh
    (boundMin,boundIdxMin,boundMax,boundIdxMax) = self.CalcBounds(mesh)

    # x component 
    #v0x= np.array(([-1,0,0]))
    #v1x= np.array(([1,0,0]))
    v0x = mesh.coordinates()[boundIdxMin[0],:]
    v1x = mesh.coordinates()[boundIdxMax[0],:]
    prob.targetsx = [tuple(v1x)]
    prob.vert_mapx = {}
    prob.vert_mapx[tuple(v0x)] = tuple(v1x)
    
    # y component 
    v0y = mesh.coordinates()[boundIdxMin[1],:]
    v1y = mesh.coordinates()[boundIdxMax[1],:]
    #v0y= np.array(([0,-1,0]))
    #v1y= np.array(([0,1,0]))
    prob.targetsy = [tuple(v1y)]
    prob.vert_mapy = {}
    prob.vert_mapy[tuple(v0y)] = tuple(v1y)
    
    # z component 
    #v0z= np.array(([0,0,-1]))
    #v1z= np.array(([0,0,1]))
    v0z = mesh.coordinates()[boundIdxMin[2],:]
    v1z = mesh.coordinates()[boundIdxMax[2],:]
    prob.targetsz= [tuple(v1z)]
    prob.vert_mapz = {}
    prob.vert_mapz[tuple(v0z)] = tuple(v1z)
  
    print v0x, " maps to " , prob.vert_mapx[ tuple(v0x)  ]
    print v0y, " maps to " , prob.vert_mapy[ tuple(v0y)  ]
    print v0z, " maps to " , prob.vert_mapz[ tuple(v0z)  ]

  def GeometryInitializations(self):
    prob = self.prob 
    mesh = prob.mesh
  
    # center
    self.CenterMesh(mesh)
  
    prob.center_coord = self.CalcMidpoint(mesh)
    #coords = mesh.coordinates()
  #
  #
    coords = mesh.coordinates()
    prob.dims = np.array([coords[:,i].max()-coords[:,i].min() for i in range(prob.nDims)])
  #  print "Center:", center_coord
  #  print "Dim:", dims
