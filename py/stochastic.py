#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import helpers
from itertools import product
import time
import ar2gas
from scipy.stats import norm

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
            
#Interpolation functions       
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
    print('Took {} seconds'.format(round((t2-t1), 2)))
    print('Inverting LHS matrix...')
    t1 = time.time()
    lhs_inv = np.linalg.inv(ones)
    t2 = time.time()
    print('Took {} seconds'.format(round((t2-t1), 2)))
    return lhs_inv

def rhs(coords, nodes, cov):
    print('Building RHS matrix...')
    krig_cov = ar2gas.compute.KrigingCovariance(1., cov)
    t1 = time.time()
    mat = []
    for pp, pg in product(coords, nodes):
        mat.append(krig_cov.compute(pg, pp))
    mat = np.array(mat).reshape((len(coords), len(nodes)))
    t2 = time.time()
    print('Took {} seconds'.format(round((t2-t1), 2)))
    return mat

def global_krig(grid, lhs_inv, rhs, var):
    t1 = time.time()
    var = np.append(var, 0)
    zeros = np.zeros((rhs.shape[0]+1, rhs.shape[1]))
    zeros[:-1,:] = rhs
    print('Computing weights...')
    t1 = time.time()
    weights = np.dot(lhs_inv, zeros)
    print('Sum of wehigths: {}'.format(round(np.sum(weights[:-1])), 2))
    t2 = time.time()
    print('Took {} seconds'.format(round((t2-t1), 2)))
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
    print('Took {} seconds'.format(round((t2-t1), 2)))
    return results

def dual_krig(grid, lhs_inv, rhs, var):
    var = np.append(var, 0)
    print('Computing weights...')
    t1 = time.time()
    weights = np.dot(lhs_inv, var)
    print('Sum of wehigths: {}'.format(round(np.sum(weights[:-1])), 2))
    t2 = time.time()
    print('Took {} seconds'.format(round((t2-t1), 2)))
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
    print('Took {} seconds'.format(round((t2-t1),2)))
    return results

def ar2gas_dual_krig(cov, x, y, z, prop, grid):
    print('Computing results by ar2gas dual kriging...')
    t1 = time.time()
    krig_cov = ar2gas.compute.KrigingCovariance(1.,cov)
    ps = ar2gas.data.PointSet(x, y, z)
    estimator = ar2gas.compute.DualKriging.OK(krig_cov, ps, prop, 0)
    tp = np.ones(grid.size_of_mask())*float('nan') if hasattr(grid, 'mask') else np.ones(grid.size())*float('nan')
    results = np.ones(grid.size())*float('nan')
    estimator.compute(grid, results, 0)
    if hasattr(grid, 'mask'):
        mask = grid.mask()
        r_idx = 0
        for idx, val in enumerate(mask):
            if val == True:
                tp[idx] = results[r_idx]
                r_idx = r_idx + 1
    else:
        tp=results
    t2 = time.time()
    print('Took {} seconds'.format(round((t2-t1),2)))
    return tp

def interpolate_variables(x, y, z, variables, codes, grid, variograms, krig_type, keep_variables, var_type, tg_prop_name, tg_grid_name):
    coords_matrix = np.vstack((x,y,z)).T
    nodes = grid.locations()
    interpolated_variables = []

    if len(variograms) == 1:
        print('Interpolating using the same covariance model for all variables')
        cov = list(variograms.values())[0]
        if krig_type != 'Ar2gas dual kriging ':
            lhs_inv = lhs(coords_matrix, cov)
            rhs_var = rhs(coords_matrix, nodes, cov)  

        for idx, v in enumerate(variables):
            rt = int(codes[idx])
            print('Interpolating RT {}'.format(rt))
            if krig_type == "Dual kriging":
                results = dual_krig(grid, lhs_inv, rhs_var, v)
                interpolated_variables.append(results)
            elif krig_type == 'Global kriging':
                results = global_krig(grid, lhs_inv, rhs_var, v)
                interpolated_variables.append(results)
            else:
                results = ar2gas_dual_krig(cov, x, y, z, v, grid)
                interpolated_variables.append(results)

            if keep_variables == '1':
                prop_name = 'interpolated_'+var_type+'_'+tg_prop_name+'_'+str(rt)
                sgems.set_property(tg_grid_name, prop_name, results.tolist())

    else:
        for idx, v in enumerate(variables):
            rt = int(codes[idx])
            print('Interpolating using one covariance model per variables')
            print('Interpolating RT {}'.format(rt))
            if krig_type != 'Ar2gas dual kriging ':
                lhs_inv = lhs(coords_matrix, variograms[rt])
                rhs_var = rhs(coords_matrix, nodes, variograms[rt])
            
            if krig_type == "Dual kriging":
                results = dual_krig(grid, lhs_inv, rhs_var, v)
                interpolated_variables.append(results)
            elif krig_type == "Global kriging":
                results = global_krig(grid, lhs_inv, rhs_var, v)
                interpolated_variables.append(results)
            else:
                results = ar2gas_dual_krig(variograms[rt], x, y, z, v, grid)
                interpolated_variables.append(results)
            
            if keep_variables == '1':
                prop_name = 'interpolated_'+var_type+'_'+tg_prop_name+'_'+str(rt)
                sgems.set_property(tg_grid_name, prop_name, results.tolist())

    print('Finished interpolating!')
    return interpolated_variables

def build_geomodel(var_type, interpolated_variables, codes, grid, tg_grid_name, tg_prop_name):
    print('Creating a geologic model...')
    
    if var_type == 'Indicators':
        
        if len(interpolated_variables) == 1:
            nn_results = helpers.nn(x, y, z, variables[0], grid)
            if keep_variables == '1':
                prop_name = 'nn_'+str(codes[0])
                sgems.set_property(tg_grid_name, prop_name, nn_results.tolist())
            proportion = sum(nn_results==1)/len(nn_results)
            print('Cutting-off interpolated indicator property in {}'.format(proportion.round(2)))
            q = np.quantile(interpolated_variables[0], (1 - proportion))
            solid = np.where(interpolated_variables[0] > q, 1, 0)
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
                    index = i.argmax(axis=0)
                    geomodel.append(float(codes[index]))
                    idx=idx+1
            sgems.set_property(tg_grid_name, tg_prop_name, geomodel)
            print('Geologic model created!')

    else:

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
    
#Simulation functions
def tbsim(mgrid, cov, nlines, nreal, seed):
    print('Simulating p-field')

    mask = mgrid.mask()
    
    tbsim = gas.compute.Tbsim.multi_realization(seed, nreal, mgrid, nlines, cov)
    results = tbsim.simulate(mgrid, nreal, 0)        
    results = [norm.cdf(lst) for lst in results]
    
    results_list = []
    for r in results:
        counter = 0
        tg_prop = np.ones(grid['nx']*grid['ny']*grid['nz']) * float('nan')
        for idx, v in enumerate(mask):
            if v == True:
                tg_prop[idx] = r[counter]
                counter = counter+1
        results_list.append(tg_prop)
                    
    return results_list
    
#block transformation functions
def sofmax_transformation(df_lst, gamma):
	
	prob_lst = np.empty(len(df_lst))

	if True in np.isnan(df_lst):
		prob_lst = np.ones(len(df_lst))*float("nan")

	else:
		exp_lst = [np.exp(-i/gamma) for i in df_lst]
		for i, exp in enumerate(exp_lst):
			prob = exp/sum(exp_lst)
			prob_lst[i] = prob

	return prob_lst

def entropy(prob_list):

	if True in np.isnan(prob_list):
		return float("nan")

	else:
		return -sum([prob*np.log2(prob) for prob in prob_list])

#################################################################################################

class stochastic: #aqui vai o nome do plugin
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
        #getting point tab variables
        tg_grid_name = self.params['gridselector']['value']
        tg_region_name = self.params['gridselector']['region']
        tg_prop_name = self.params['lineEdit']['value']
        pt_grid_name = self.params['gridselectorbasic']['value']
        pt_props_name = self.params['orderedpropertyselector']['value'].split(';')
        keep_variables = self.params['checkBox_2']['value']
        krig_type = self.params['comboBox_2']['value']
        fx, fy, fz = int(self.params['fx']['value']), int(self.params['fy']['value']), int(self.params['fz']['value'])
        codes = [v.split('_')[-1] for v in pt_props_name]
        
        #getting block properties
        re_use = self.params['checkBox_3']['value']
        gamma = float(self.params['doubleSpinBox']['value'])
        if re_use == '1 ':
            grid_grid_name = self.params['propertyselector']['grid']
            grid_entropy = self.params['propertyselector']['property']
            grid_prob_name = self.params['gridselectorbasic_2']['value']
            grid_pobs_props_names = self.params['orderedpropertyselector_2']['value'].split(';')
        
        #p-field tab variables
        n_reals = self.params['spinBox']['value'] 
        seed = self.params['spinBox_3']['value'] 
        nlines = self.params['spinBox_2']['value'] 
        
        #getting variograms
        p = self.params
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
            varg_lst = helpers.ar2gemsvarwidget_to_ar2gascovariance(p)
            if len(varg_lst) == 1:
                varg_lst = varg_lst * len(codes)
                variograms = {}
                variograms[0] = varg_lst[0]
            else:
                variograms = dict(zip(codes, varg_lst))
                
        pfieldvar = helpers.singlear2gemsvarwidget_to_ar2gascovariance(p)
        
        #interpolating to all grid nodes if not re-using properties
        if re_use == '0':
            print('Not re-using properties')
            a2g_grid = helpers.ar2gemsgrid_to_ar2gasgrid(tg_grid_name, tg_region_name)
            
            variables = []
            nan_filters = []
            for v in pt_props_name:
                values = np.array(sgems.get_property(pt_grid_name, v))
                nan_filter = np.isfinite(values)
                filtered_values = values[nan_filter]
                nan_filters.append(nan_filter)
                variables.append(filtered_values)
            nan_filter = np.product(nan_filters, axis=0)
            nan_filter = nan_filter == 1

            x, y, z = np.array(sgems.get_X(pt_grid_name))[nan_filter], np.array(sgems.get_Y(pt_grid_name))[nan_filter], np.array(sgems.get_Z(pt_grid_name))[nan_filter]
            
            var_type = ''
            interpolated_variables = interpolate_variables(x, y, z, variables, codes, a2g_grid, variograms, krig_type, keep_variables, var_type, tg_prop_name, tg_grid_name)
            interpolated_variables = np.array(interpolated_variables)
            
            #calculating probs
            print('Calculating probabilities...')
            probs_matrix = np.array([sofmax_transformation(sds, gamma) for sds in interpolated_variables.T])
            probs_matrix = probs_matrix.T
            
            if keep_variables == '1':
                for i, p in enumerate(probs_matrix):
                   sgems.set_property(tg_grid_name, pt_props_name[i]+'_gamma_'+str(gamma), p.tolist())
               
            #calculating entropy
            print('Calculating entropy...')
            entropies = [entropy(probs) for probs in probs_matrix.T]
            entropies = (entropies - min(entropies))/(max(entropies) - min(entropies))
            if keep_variables == '1':
                sgems.set_property(tg_grid_name, 'entropy_gamma_'+str(gamma), entropies.tolist())
            
            print('Done!')
          
        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "stochastic" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["stochastic"] #aqui vai o nome do plugin