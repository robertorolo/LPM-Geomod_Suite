#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import helpers
from scipy.interpolate import NearestNDInterpolator
import ar2gas
from itertools import product


#################################################################################################

#defina aqui as funcoes que serao utilizadas no codigo

#exibe os valores do dicionario de parametros
def read_params(a,j=''):
    if j=='':
        print("### Printing GUI parameters ###")
    for i in a:
        if (type(a[i])!=type({'a':1})):
            print(j+"['"+str(i)+"']="+str(a[i])+" type: "+str(type(a[i])))
        else:
            read_params(a[i],j+"['"+str(i)+"']")

def marching_cubes(grid):
    nx, ny, nz = grid.dim()[0], grid.dim()[1], grid.dim()[2]

    range_x, range_y, range_z = [l for l in range(0,nx-1)], [l for l in range(0,ny-1)], [l for l in range(0,nz-1)]

    indices_list = []
    
    for k, j, i in product(range_z, range_y, range_x):
    
        cube = np.array([[i,j,k],[i,j+1,k],[i+1,j+1,k],[i+1,j,k],[i,j,k+1],[i,j+1,k+1],[i+1,j+1,k+1],[i+1,j,k+1]])
        indices = [helpers.ijk_in_n(grid, e[0], e[1], e[2]) for e in cube]
        indices_list.append(indices)
    
    return indices_list

def refinement_zone(grid, geomodel):
    refinement_prop = geomodel.copy()
    indices_list = marching_cubes(grid)
    for indice_list in indices_list:
        cats = geomodel[indice_list]
        if np.unique(cats).size is not 1:
            refinement_prop[indice_list] = -999

    return refinement_prop, indices_list

def nn(x, y, z, var, grid):
    nan_mask = np.isfinite(var)
    points_array = np.vstack((x,y,z)).T
    knn = NearestNDInterpolator(points_array[nan_mask], var[nan_mask])
    grids_points_array = grid.locations()
    results = knn(grids_points_array)

    return results

def dual_kriging(cov, x, y, z, prop, grid):
    krig_cov = ar2gas.compute.KrigingCovariance(1.,cov)
    ps = ar2gas.data.PointSet(np.array(x), np.array(y), np.array(z))
    estimator = ar2gas.compute.DualKriging.OK(krig_cov, ps, np.array(prop), 0)
    target_prop = np.ones(grid.size())*float('nan')
    estimator.compute(grid, target_prop, 0)

    return target_prop

#################################################################################################

class deterministic: #aqui vai o nome do plugin
    def __init__(self):
        pass

#################################################################################################

    def initialize(self, params):
        self.params = params
        
        #imprimindo o dicionario de parametros
        print("dicionario de parametros: ", params)

        #executando a funcao exibe os valores do dicionario de parametros
        read_params(params) #para nao printar comente essa linha

        return True

#################################################################################################

    def execute(self):
      
        #aqui vai o codigo

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "deterministic" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["deterministic"] #aqui vai o nome do plugin