#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import helpers
import ar2gas
import time
from scipy.spatial.distance import cdist
import math

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
    
    if nz != 1:
        range_x, range_y, range_z = [l for l in range(0,nx-1)], [l for l in range(0,ny-1)], [l for l in range(0,nz-1)]
        indices_list = []       
        for k, j, i in product(range_z, range_y, range_x):       
            cube = np.array([[i,j,k],[i,j+1,k],[i+1,j+1,k],[i+1,j,k],[i,j,k+1],[i,j+1,k+1],[i+1,j+1,k+1],[i+1,j,k+1]])
            indices = [helpers.ijk_in_n(grid, e[0], e[1], e[2]) for e in cube]
            indices_list.append(indices)

    else:
        range_x, range_y = [l for l in range(0,nx-1)], [l for l in range(0,ny-1)]
        indices_list = []        
        for j, i in product(range_y, range_x):       
            cube = np.array([[i,j,0],[i,j+1,0],[i+1,j+1,0],[i+1,j,0],[i,j,0],[i,j+1,0],[i+1,j+1,0],[i+1,j,0]])
            indices = [helpers.ijk_in_n(grid, e[0], e[1], e[2]) for e in cube]
            indices_list.append(indices)
    
    return indices_list

def refinement_zone(grid, geomodel):
    geomodel = np.array(geomodel)
    refinement_prop = geomodel.copy()
    indices_list = marching_cubes(grid)
    for indice_list in indices_list:
        cats = geomodel[indice_list]
        if np.unique(cats).size is not 1:
            refinement_prop[indice_list] = -999

    return refinement_prop, indices_list

def build_geomodel(interpolated_variables, codes, grid, tg_grid_name, tg_prop_name):
    print('Creating a geologic model...')

    if len(interpolated_variables) == 1:
        solid = np.where(interpolated_variables[0] < 0, 1, 0)
        sgems.set_property(tg_grid_name, 'rt_{}'.format(codes[0]), solid.tolist())

    else:
        int_variables = np.array(interpolated_variables).T
        geomodel = []
        idx = 0
        for i in int_variables:
            if np.isnan(i).all():
                geomodel.append(float('nan'))
                idx=idx+1
            else:
                index = i.argmin(axis=0)
                geomodel.append(float(codes[index]))
                idx=idx+1
        sgems.set_property(tg_grid_name, tg_prop_name, geomodel)
        print('Geologic model created!')
            
    return geomodel

def build_refined_geomodel(interpolated_variables, downscaled_props, codes, var_type):
    print('Creating a geologic model...')
    final_geomodel = []
    idx = 0
    for i in np.array(interpolated_variables).T:
        if np.isnan(i).all():
            final_geomodel.append(float(downscaled_props[1][idx]))
            idx=idx+1
        else:
            index = i.argmin(axis=0) if var_type == 'Signed Distances' else i.argmax(axis=0)
            final_geomodel.append(float(codes[index]))
            idx=idx+1

    return final_geomodel

def kernel(function, nugget, support, dist):
    support = 1.0/support
    er=support*dist
    if function is 'gaussian':
        return (1.0 - nugget) * np.exp(-er**1.8) if dist is not 0 else nugget
    if function is 'spherical':
        return (1.0 - nugget) * (1.0 - ( 1.5 * er - 0.5 * (er) ** 3)) if dist is not 0 else nugget

#################################################################################################

class data_conditioning_uncertainty: #aqui vai o nome do plugin
    def __init__(self):
        pass

#################################################################################################

    def initialize(self, params):
        self.params = params
        
        #imprimindo o dicionario de parametros
        #print("dicionario de parametros: ", params)

        #executando a funcao exibe os valores do dicionario de parametros
        #read_params(params) #para nao printar comente essa linha

        return True

#################################################################################################

    def execute(self):
      
        #aqui vai o codigo
        #getting variables

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "data_conditioning_uncertainty" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["data_conditioning_uncertainty"] #aqui vai o nome do plugin
