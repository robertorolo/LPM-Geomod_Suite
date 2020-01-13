#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import helpers
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial.distance import cdist
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

def lhs(coords_matrix, cov):
    print('Calculating LHS matrix...')
    krig_cov = ar2gas.compute.KrigingCovariance(1., cov)
    n = len(coords_matrix)
    ones = np.ones((n+1, n+1))
    lhs = krig_cov.lhs(coords_matrix)
    print(lhs)
    ones[:-1,:-1] = lhs
    ones[n,n] = 0
    print('Inverting LHS matrix...')
    lhs_inv = np.linalg.inv(ones)
    return lhs_inv

def global_krig(coords, grid, cov, lhs_inv, var):
    krig_cov = ar2gas.compute.KrigingCovariance(1., cov)
    var = np.append(var, 0)
    if hasattr(grid, 'mask'):
        mask = grid.mask()
        mask = np.where(mask == True, 1, 0)
    else:
        mask = np.ones(grid.size())
    nodes = grid.locations()
    result = np.ones(grid.size())*float('nan')
    for idx, maskval in enumerate(mask):
        if idx%1000 == 0:
            print('Interpolating node {}'.format(idx))
        if maskval == 1:
            rhs = krig_cov.rhs(coords, nodes[idx])
            rhs = np.append(rhs, 0)
            w = np.dot(rhs, lhs_inv).T
            z = np.dot(w, var)
            result[idx] = z
    return result

#################################################################################################

class deterministic: #aqui vai o nome do plugin
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
        var_type = self.params['comboBox']['value']
        tg_grid_name = self.params['gridselector']['value']
        tg_region_name = self.params['gridselector']['region']
        tg_prop_name = self.params['lineEdit']['value']
        keep_variables = self.params['checkBox_2']['value']
        props_grid_name = self.params['gridselectorbasic']['value']
        n_var = self.params['orderedpropertyselector']['count']
        var_names = self.params['orderedpropertyselector']['value'].split(';')
        codes = [v.split('_')[-1] for v in var_names]

        #getting variograms
        print('Getting variogram models')
        use_model_file = self.params['checkBox']['value'] 
        if use_model_file == '1':
            path = self.params['filechooser']['value']
            variograms = helpers.modelfile_to_ar2gasmodel(path)
            if len(variograms) == 1:
                values_covs = list(variograms.values())
                varg_lst = values_covs * len(codes)
                variograms = {}
                variograms[0] = varg_lst[0]
        else:
            p = self.params
            varg_lst = helpers.ar2gemsvarwidget_to_ar2gascovariance(p)
            if len(varg_lst) == 1:
                varg_lst = varg_lst * len(codes)
                variograms = {}
                variograms[0] = varg_lst[0]
            else:
                variograms = dict(zip(codes, varg_lst))
            
        #interpolating variables
        a2g_grid = helpers.ar2gemsgrid_to_ar2gasgrid(tg_grid_name, tg_region_name)
        
        variables = []
        nan_filters = []
        for v in var_names:
            values = np.array(sgems.get_property(props_grid_name, v))
            nan_filter = np.isfinite(values)
            values = np.array(sgems.get_property(props_grid_name, v))[nan_filter]
            nan_filters.append(nan_filter)
            variables.append(values)
        nan_filter = np.product(nan_filters, axis=0)
        nan_filter = nan_filter == 1
        print(nan_filter)

        x, y, z = np.array(sgems.get_X(props_grid_name))[nan_filter], np.array(sgems.get_Y(props_grid_name))[nan_filter], np.array(sgems.get_Z(props_grid_name))[nan_filter]
        coords_matrix = np.vstack((x,y,z)).T  

        if len(variograms) == 1:
            print('Interpolating using the same covarinace model for all variables')
            cov = list(variograms.values())[0]
            lhs_inv = lhs(coords_matrix, cov)

            for idx, v in enumerate(var_names):
                rt = codes[idx]
                print('Interpolating RT {}'.format(rt))
                results = global_krig(coords_matrix, a2g_grid, cov, lhs_inv, variables[idx])
                if keep_variables == '1':
                    prop_name = 'interpolated_'+var_type+'_'+tg_prop_name+'_'+str(rt)
                    sgems.set_property(tg_grid_name, prop_name, results.tolist())

        else:
            for idx, v in enumerate(var_names):
                rt = codes[idx]
                print('Interpolating using one covarinace model per variables')
                print('Interpolating RT {}'.format(rt))
                lhs_inv = lhs(coords_matrix, variograms[idx])
                results = global_krig(coords_matrix, a2g_grid, variograms[rt], lhs_inv, variables[idx])
                if keep_variables == '1':
                    prop_name = 'interpolated_'+var_type+'_'+tg_prop_name+'_'+str(rt)
                    sgems.set_property(tg_grid_name, prop_name, results.tolist())

        print('Finished!')

        #creating a geologic model


        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "deterministic" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["deterministic"] #aqui vai o nome do plugin