#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
#np.set_printoptions(threshold=np.inf)
import helpers
import ar2gas
from itertools import product
import time
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import norm
from sklearn.metrics import confusion_matrix

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

def ar2gas_dual_krig(cov, x, y, z, prop, grid):
    
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

def perturb_sd(x, y, z, variables, codes, sim_variograms, seed, n_lines, p_factor):
    pset = ar2gas.data.PointSet(x, y, z)
    bb = pset.bounding_box()
    factors = []
    for idx, v in enumerate(variables):
        rt = int(codes[idx])
        if len(sim_variograms) == 1:
            cov = list(sim_variograms.values())[0] 
        else:
            cov = sim_variograms[rt]
        tbsim = ar2gas.compute.Tbsim(bb, n_lines, cov, seed, 0)
        results = np.ones(pset.size()) * float('nan')
        tbsim.simulate(pset, results)
        inverse = norm.cdf(results)
        shape = inverse.shape
        inverse = inverse.reshape(-1, 1)
        scaler = MinMaxScaler(feature_range=(1-p_factor, 1+p_factor))
        scaled = scaler.fit_transform(inverse).reshape(shape)
        factors.append(scaled)

    return np.array(variables) * np.array(factors)

def interpolate_variables(x, y, z, variables, codes, grid, vargs, keep_variables, var_type, tg_prop_name, tg_grid_name):
    print('Computing results by ar2gas dual kriging...')

    interpolated_variables = []

    if len(vargs) == 0:
        print('Interpolating using the same covariance model for all variables')

        for idx, v in enumerate(variables):
            rt = int(codes[idx])
            print('Interpolating RT {}'.format(rt))
            if var_type == 'Indicators':
                mask = v==1
            else:
                mask = v < 0
            var_range = helpers.max_dist(x, y, z, mask, grid)
            print('Estimated range is {}'.format(var_range))
            nugget = 0.01
            cov = [ar2gas.compute.Covariance.nugget(nugget), ar2gas.compute.Covariance.gaussian(1-nugget, var_range, var_range, var_range, 0., 0., 0.)]
            vargs[rt] = cov
            results = ar2gas_dual_krig(cov, x, y, z, v, grid)
            interpolated_variables.append(results)

            if keep_variables == '1':
                prop_name = 'interpolated_'+var_type+'_'+tg_prop_name+'_'+str(rt)
                sgems.set_property(tg_grid_name, prop_name, results.tolist())

    elif len(vargs) == 1:
        print('Interpolating using the same covariance model for all variables')
        cov = list(vargs.values())[0] 

        for idx, v in enumerate(variables):
            rt = int(codes[idx])
            print('Interpolating RT {}'.format(rt))
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
            results = ar2gas_dual_krig(vargs[rt], x, y, z, v, grid)
            interpolated_variables.append(results)
            
            if keep_variables == '1':
                prop_name = 'interpolated_'+var_type+'_'+tg_prop_name+'_'+str(rt)
                sgems.set_property(tg_grid_name, prop_name, results.tolist())

    print('Finished interpolating!')
    return interpolated_variables

def build_geomodel(var_type, interpolated_variables, codes, grid, tg_grid_name, tg_prop_name, ids, acceptance, rt_prop, num_models):
    print('Creating a geologic model...')
    
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
    cm = confusion_matrix(rt_prop, np.array(geomodel)[ids], normalize='true')
    diag = np.diagonal(cm)
    m = np.mean(diag)
    if m < acceptance:
        pass
    else:
        sgems.set_property(tg_grid_name, tg_prop_name, geomodel)
        num_models.append('x')
            
    return geomodel

def sd_to_cat(variables, codes):
    target = np.zeros(len(variables[0]))
    for idx, v in enumerate(variables):
        target[v < 0] = codes[idx]
    return target

#################################################################################################

class sd_perturbation_uncertainty: #aqui vai o nome do plugin
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
        tg_grid_name = self.params['gridselector']['value']
        tg_region_name = self.params['gridselector']['region']
        tg_prop_name = self.params['lineEdit']['value']
        keep_variables = self.params['checkBox_2']['value']
        props_grid_name = self.params['gridselectorbasic']['value']
        var_names = self.params['orderedpropertyselector']['value'].split(';')
        codes = [v.split('_')[-1] for v in var_names]
        var_type = 'Signed Distances'

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

        #sim variables
        n_reals = int(self.params['spinBox']['value'])
        n_lines = int(self.params['spinBox_2']['value'])
        p_factor = float(self.params['doubleSpinBox']['value'])
        seed = int(self.params['spinBox_3']['value'])
        acceptance = float(self.params['doubleSpinBox_2']['value'])

        #getting variograms for interpolation
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
            if len(varg_lst) == 0:
                variograms = {}
            if len(varg_lst) == 1:
                varg_lst = varg_lst * len(codes)
                variograms = {}
                variograms[0] = varg_lst[0]

            else:
                variograms = dict(zip(codes, varg_lst))

        #getting variograms for simulation
        use_model_file = self.params['checkBox_3']['value']
        if use_model_file == '1':
            path = self.params['filechooser_2']['value']
            sim_variograms = helpers.modelfile_to_ar2gasmodel(path)
            if len(variograms) == 1:
                values_covs = list(sim_variograms.values())
                varg_lst = values_covs * len(codes)
                sim_variograms = {}
                sim_variograms[0] = varg_lst[0]
        else:
            p = self.params
            varg_lst = helpers.ar2gemsvarwidget_to_ar2gascovariance_sim(p)
            if len(varg_lst) == 0:
                sim_variograms = {}
            if len(varg_lst) == 1:
                varg_lst = varg_lst * len(codes)
                sim_variograms = {}
                sim_variograms[0] = varg_lst[0]

            else:
                sim_variograms = dict(zip(codes, varg_lst))
            
        #interpolating variables
        num_models = []
        for idx in range(n_reals):
            print('Working on realization {}...'.format(idx+1))

            tg_prop_name_temp = tg_prop_name+'_real_{}'.format(idx)

            perturbed_variables = perturb_sd(x, y, z, variables, codes, sim_variograms, seed, n_lines, p_factor)
            seed = seed + 1
            
            a2g_grid = helpers.ar2gemsgrid_to_ar2gasgrid(tg_grid_name, tg_region_name)
            
            interpolated_variables = interpolate_variables(x, y, z, perturbed_variables, codes, a2g_grid, variograms, keep_variables, var_type, tg_prop_name_temp, tg_grid_name)

            #creating a geologic model
            ids = [sgems.get_closest_nodeid(tg_grid_name, xi, yi, zi) for xi, yi, zi in zip(x,y,z)]
            rt_prop = sd_to_cat(variables, codes)
            build_geomodel(var_type, interpolated_variables, codes, a2g_grid, tg_grid_name, tg_prop_name_temp, ids, acceptance, rt_prop, num_models)

        print('{} geologic models accepted!'.format(len(num_models)))
        print('Finished!')

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "sd_perturbation_uncertainty" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["sd_perturbation_uncertainty"] #aqui vai o nome do plugin
