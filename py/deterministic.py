#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
np.set_printoptions(threshold=np.inf)
import helpers
from scipy.interpolate import NearestNDInterpolator
import ar2gas
from itertools import product
import time

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
            indices = [ijk_in_n(grid, e[0], e[1], e[2]) for e in cube]
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

def lhs(coords_matrix, cov, global_krig=False):
    print('Calculating LHS matrix...')
    t1 = time.time()
    krig_cov = ar2gas.compute.KrigingCovariance(1., cov)
    lhs = krig_cov.lhs(coords_matrix)
    n = len(coords_matrix)
    ones = np.ones((n+1, n+1))
    ones[:-1,:-1] = lhs
    ones[n,n] = 0
    t2 = time.time()
    print('Took {} seconds'.format(t2-t1))
    print('Inverting LHS matrix...')
    t1 = time.time()
    lhs_inv = np.linalg.inv(ones)
    t2 = time.time()
    print('Took {} seconds'.format(t2-t1))
    return lhs_inv

def matrices_operations(krig_cov, coords, node): #can be used in paralelization
    rhs = krig_cov.rhs(coords, node) #this method is not working
    rhs = np.append(rhs, 0)
    w = np.dot(rhs, lhs_inv).T
    z = np.dot(w, var)
    return z

def rhs(coords, nodes, cov):
    print('Building RHS matrix...')
    krig_cov = ar2gas.compute.KrigingCovariance(1., cov)
    t1 = time.time()
    n = len(nodes)
    i = 0
    mat = []
    for pp, pg in product(coords, nodes):
        mat.append(krig_cov.compute(pg, pp))
    mat = np.array(mat).reshape((len(coords), len(nodes)))
    t2 = time.time()
    print('Took {} seconds'.format(t2-t1))
    return mat

def global_krig(grid, lhs_inv, rhs, var):
    t1 = time.time()
    var = np.append(var, 0)
    zeros = np.zeros((rhs.shape[0]+1, rhs.shape[1]))
    zeros[:-1,:] = rhs
    print('Computing weights...')
    t1 = time.time()
    weights = np.dot(lhs_inv, zeros)
    print('Sum of wehigths: {}'.format(np.sum(weights[:-1])))
    t2 = time.time()
    print('Took {} seconds'.format(t2-t1))
    print('Computing results by global kriging...')
    t1 = time.time()
    partial_results = np.dot(var, weights)

    if hasattr(grid, 'mask'):
        mask = grid.mask()

        results = np.ones(len(mask))*float('nan')
        mg_idx = 0
        for idx, maskval in enumerate(mask):
            if maskval == True:
                results[idx] = partial_results[mg_idx]
                mg_idx = mg_idx+1
    else:
        results = partial_results
    
    t2 = time.time()
    print('Took {} seconds'.format(t2-t1))
    return results

def dual_krig(grid, lhs_inv, rhs, var):
    var = np.append(var, 0)
    print('Computing weights...')
    t1 = time.time()
    weights = np.dot(lhs_inv, var)
    print('Sum of wehigths: {}'.format(np.sum(weights[:-1])))
    t2 = time.time()
    print('Took {} seconds'.format(t2-t1))
    print('Computing results by dual kriging...')
    t1 = time.time()
    partial_results = np.dot(rhs.T, weights[:-1]) + weights[-1:]

    if hasattr(grid, 'mask'):
        mask = grid.mask()

        results = np.ones(len(mask))*float('nan')
        mg_idx = 0
        for idx, maskval in enumerate(mask):
            if maskval == True:
                results[idx] = partial_results[mg_idx]
                mg_idx = mg_idx+1
    else:
        results = partial_results
    
    t2 = time.time()
    print('Took {} seconds'.format(t2-t1))
    return results

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
        krig_type = self.params['comboBox_2']['value']

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
        nodes = a2g_grid.locations()
        
        variables = []
        nan_filters = []
        for v in var_names:
            values = np.array(sgems.get_property(props_grid_name, v))
            nan_filter = np.isfinite(values)
            filtered_values = values[nan_filter]
            nan_filters.append(nan_filter)
            variables.append(filtered_values)
        nan_filter = np.product(nan_filters, axis=0)
        nan_filter = nan_filter == 1

        x, y, z = np.array(sgems.get_X(props_grid_name))[nan_filter], np.array(sgems.get_Y(props_grid_name))[nan_filter], np.array(sgems.get_Z(props_grid_name))[nan_filter]
        coords_matrix = np.vstack((x,y,z)).T

        if len(variograms) == 1:
            print('Interpolating using the same covarinace model for all variables')
            cov = list(variograms.values())[0]
            lhs_inv = lhs(coords_matrix, cov)
            rhs_var = rhs(coords_matrix, nodes, cov)

            for idx, v in enumerate(var_names):
                rt = codes[idx]
                print('Interpolating RT {}'.format(rt))
                if krig_type == "Dual kriging":
                    results = dual_krig(a2g_grid, lhs_inv, rhs_var, variables[idx])
                else:
                    results = global_krig(a2g_grid, lhs_inv, rhs_var, variables[idx])
                if keep_variables == '1':
                    prop_name = 'interpolated_'+var_type+'_'+tg_prop_name+'_'+str(rt)
                    sgems.set_property(tg_grid_name, prop_name, results.tolist())

        else:
            for idx, v in enumerate(var_names):
                rt = codes[idx]
                print('Interpolating using one covarinace model per variables')
                print('Interpolating RT {}'.format(rt))
                lhs_inv = lhs(coords_matrix, variograms[rt])
                rhs_var = rhs(coords_matrix, nodes, variograms[rt])
                if krig_type == "Dual kriging":
                    results = dual_krig(a2g_grid, lhs_inv, rhs_var, variables[idx])
                else:
                    results = global_krig(a2g_grid, lhs_inv, rhs_var, variables[idx])
                if keep_variables == '1':
                    prop_name = 'interpolated_'+var_type+'_'+tg_prop_name+'_'+str(rt)
                    sgems.set_property(tg_grid_name, prop_name, results.tolist())

        print('Finished interpolating!')

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