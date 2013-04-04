"""
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
"""
# Defines boundary conditions, etc for microdomain
# See notes from 120514_grids.tex
from dolfin import *
from params import *
from Domain import *

def boundary(x,on_boundary):
  return on_boundary


class CellularUnitDomain_TnC(Domain):
  def __init__(self,type):
    super(CellularUnitDomain_TnC,self).__init__(type)

    problem = self.problem
    problem.name = "Cellular"

  def Setup(self):
    # mesh
    problem = self.problem
    problem.mesh = UnitCube(8,8,8)

    # resizing mesh to be 0.05x0.05x0.05 [um]
    problem.mesh.coordinates()[:] *= 0.05/2. # since centered at origin


    if(self.type=="scalar"):
        problem.V = FunctionSpace(problem.mesh,"CG",1)
    elif(self.type=="field"):
        problem.V = VectorFunctionSpace(problem.mesh,"CG",1)

    # geom
    self.CalcGeom(problem)

  # bcs
  def AssignBC(self):
    problem = self.problem

    if(self.type=="scalar"):
        u0 = Constant(0.)
        #u1 = Expression("1 + x[0]*x[0] + x[1]*x[1]")
        # I don't think I want any specific BCs here for homog prob.
        u1 = Constant(1.)
    elif(self.type=="field"):
        u0 = Constant((0.,0,0.))
        u1 = Constant((1.,1.,1.))


    print "Not sure if I should be using any Dirichlet BC"
    bc1 = DirichletBC(problem.V,u1,boundary)
    # neum

    problem.bcs = [bc1]


