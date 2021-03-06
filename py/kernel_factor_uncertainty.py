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
import itertools
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import MinMaxScaler

#################################################################################################

#defina aqui as funcoes que serao utilizadas no codigo

#################################################################################################
#RBF functions and classes

def coordinates_transform(coords, major_med, major_min, azimuth, dip, rake):

    if azimuth >= 0 and azimuth <=270:
        alpha = math.radians(90-azimuth)
    else:
        alpha = math.radians(450-azimuth)
    beta = -math.radians(dip)
    phi = math.radians(rake)
    
    rot_matrix = np.zeros((3,3))

    rot_matrix[0,0] = math.cos(beta)*math.cos(alpha)
    rot_matrix[0,1] = math.cos(beta)*math.sin(alpha)
    rot_matrix[0,2] = -math.sin(beta)
    rot_matrix[1,0] = major_med*(-math.cos(phi)*math.sin(alpha)+math.sin(phi)*math.sin(beta)*math.cos(alpha))
    rot_matrix[1,1] = major_med*(math.cos(phi)*math.cos(alpha)+math.sin(phi)*math.sin(beta)*math.sin(alpha))
    rot_matrix[1,2] = major_med*(math.sin(phi)*math.cos(beta))
    #rot_matrix[2,0] = major_min*(math.sin(phi)*math.sin(alpha)+math.cos(phi)*math.sin(beta)*math.cos(alpha))
    rot_matrix[2,0] = major_min*(math.sin(phi)*math.sin(alpha)+math.cos(phi)*math.cos(beta)*math.cos(alpha))
    rot_matrix[2,1] = major_min*(-math.sin(phi)*math.cos(alpha)+math.cos(phi)*math.sin(beta)*math.sin(alpha))
    rot_matrix[2,2] = major_min*(math.cos(phi)*math.cos(beta))

    return np.array([np.dot(rot_matrix, i) for i in coords])

def kernel(function, nugget, support, dist):
    support = 1.0/support
    er=support*dist
    if function == 'gaussian':
        return (1.0 - nugget) * np.exp(-er**1.8) if dist is not 0 else nugget
    if function == 'spherical':
        return (1.0 - nugget) * (1.0 - ( 1.5 * er - 0.5 * (er) ** 3)) if dist is not 0 else nugget

class RBF:

    def __init__(self, x, y, z, var, function='gaussian', nugget=0.01, support=0, major_med=1, major_min=1, azimuth=0, dip=0, rake=0):
        self.x = x
        self.y = y
        self.z = z
        self.var = var
        self.function = function
        self.nugget = nugget        
        self.support = support
        self.major_med = major_med
        self.major_min = major_min
        self.azimuth = 0
        self.dip = 0
        self.rake = rake
        self.weights = None

        X = np.array([x, y, z]).T
        self.X = coordinates_transform(X, major_med, major_min, azimuth, dip, rake)

    def train(self, dcf):
        print('Training weights...')
        self.supports = self.support * dcf
        self.supports = self.supports.reshape(len(self.supports), 1)
        print('max support: ', np.max(self.supports), 'min support: ', np.min(self.supports))
        dist_mat = cdist(self.X, self.X)
        int_mat = kernel(self.function, self.nugget, self.supports, dist_mat)
        if np.linalg.det(int_mat) != 0:
            weights = np.dot(np.linalg.inv(int_mat), self.var)
        else:
            print('Det equal zero. Calculating pseudo inverse...')
            weights = np.dot(np.linalg.pinv(int_mat), self.var)
        self.weights = weights

    def predict(self, a2ggrid):
        print('Predicting...')

        x = a2ggrid.locations()
        x = coordinates_transform(x, self.major_med, self.major_min, self.azimuth, self.dip, self.rake)

        dist_mat = cdist(self.X, x)
        int_mat = kernel(self.function, self.nugget, self.supports, dist_mat)
        results = np.dot(int_mat.T, self.weights)

        #accounting for the mask
        tp = np.ones(a2ggrid.size_of_mask())*float('nan') if hasattr(a2ggrid, 'mask') else np.ones(a2ggrid.size())*float('nan')
        if hasattr(a2ggrid, 'mask'):
            mask = a2ggrid.mask()
            r_idx = 0
            for idx, val in enumerate(mask):
                if val == True:
                    tp[idx] = results[r_idx]
                    r_idx = r_idx + 1
        else:
            tp=results
            
        return tp

#################################################################################################

#exibe os valores do dicionario de parametros
def read_params(a,j=''):
    if j=='':
        print("### Printing GUI parameters ###")
    for i in a:
        if (type(a[i])!=type({'a':1})):
            print(j+"['"+str(i)+"']="+str(a[i])+" type: "+str(type(a[i])))
        else:
            read_params(a[i],j+"['"+str(i)+"']")

def sqrt_coef(vertice, p):
    h = vertice[0]
    k = vertice[1]
    a = (p[1]-k)/np.sqrt(p[0]-h)
    
    return a, h, k

def line_coef(p1, p2):
    x1, y1 = p1[0], p1[1]
    x2, y2 = p2[0], p2[1]
    a = (y2 - y1) / (x2 - x1)
    b = y1 - a * x1     
    
    return a, b

def dc_parametrization(sd, function, f_min_o, f_min_p):
    param = {}
    sd_min = np.min(sd)
    sd_max = np.max(sd)
    if function == 'linear':
        ap, bp = line_coef([sd_min, f_min_p], [0, 1])
        ao, bo = line_coef([sd_max, f_min_o], [0, 1])
        param['pessimistic'] = ap*sd+bp
        param['optimistic'] = ao*sd+bo
        param['intermediate'] = np.ones(len(sd))
    else:
        x_neg = sd[sd<0]
        x_pos = sd[sd>0]
        
        a_pos_p, h_pos_p, k_pos_p = sqrt_coef(vertice=[0, 1], p=[sd_max, f_min_p])
        a_neg_p, h_neg_p, k_neg_p = sqrt_coef(vertice=[0, 1], p=[np.abs(sd_min), f_min_p])

        a_pos_o, h_pos_o, k_pos_o = sqrt_coef(vertice=[0, 1], p=[sd_max, f_min_o])
        a_neg_o, h_neg_o, k_neg_o = sqrt_coef(vertice=[0, 1], p=[np.abs(sd_min), f_min_o])
        
        y_pos = a_pos_o * np.sqrt(x_pos-h_pos_o) + k_pos_o
        y_neg = -a_neg_o * np.sqrt(np.abs(x_neg-h_neg_o)) + k_neg_o
        yo = np.concatenate([y_neg, y_pos])

        y_pos = -a_pos_p * np.sqrt(x_pos-h_pos_p) + k_pos_p
        y_neg = a_neg_p * np.sqrt(np.abs(x_neg-h_neg_p)) + k_neg_p
        yp = np.concatenate([y_neg, y_pos])

        param['pessimistic'] = yp
        param['optimistic'] = yo
        param['intermediate'] = np.ones(len(sd))
    
    return param

def combinate(codes):
    tps = ['pessimistic', 'intermediate', 'optimistic']
    return list(itertools.product(tps, repeat=len(codes)))

def build_geomodels(interpolated_variables, codes, tg_grid_name, tg_prop_name, ids, acceptance, rt_prop):
    #implicit binary modeling
    if len(interpolated_variables) == 1:
        for t in interpolated_variables[int(codes[0])]:
            solid = np.where(interpolated_variables[int(codes[0])][t] < 0, 1, 0)
            sgems.set_property(tg_grid_name, 'rt_{}_{}'.format(codes[0], t), solid.tolist())

    #multicategorical
    else:
        accepted_combs = []
        num_models = 0
        combinations = combinate(codes)
        for sc, c in enumerate(combinations):
            int_variables = []
            for idx, r in enumerate(codes):
                int_variables.append(interpolated_variables[int(r)][c[idx]])
            
            geomodel = []
            idx = 0
            for i in np.array(int_variables).T:
                if np.isnan(i).all():
                    geomodel.append(float('nan'))
                    idx=idx+1
                else:
                    index = i.argmin(axis=0)
                    geomodel.append(float(codes[index]))
                    idx=idx+1
            
            cm = confusion_matrix(rt_prop, np.array(geomodel)[ids], normalize='true')
            diag = np.diagonal(cm)
            print(diag)
            m = np.mean(diag)
            #if m < acceptance:
            if np.any(diag < acceptance):
                pass
            else:
                sgems.set_property(tg_grid_name, tg_prop_name+'_'+str(sc), geomodel)
                num_models = num_models + 1
                accepted_combs.append(c)
        print('{} geologic models accepted!'.format(num_models))
        print(accepted_combs)

def sd_to_cat(variables, codes):
    target = np.zeros(len(variables[0]))
    for idx, v in enumerate(variables):
        target[v < 0] = codes[idx]
    return target

def mask_var(a2ggrid, results):
    tp = np.ones(a2ggrid.size_of_mask())*float('nan') if hasattr(a2ggrid, 'mask') else np.ones(a2ggrid.size())*float('nan')
    if hasattr(a2ggrid, 'mask'):
        mask = a2ggrid.mask()
        r_idx = 0
        for idx, val in enumerate(mask):
            if val == True:
                tp[idx] = results[r_idx]
                r_idx = r_idx + 1
    else:
        tp=results
        
    return tp

#################################################################################################

class kernel_factor_uncertainty: #aqui vai o nome do plugin
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
        
        dcf_param = self.params['comboBox_2']['value']
        f_min_o = float(self.params['doubleSpinBox']['value'])
        f_min_p = float(self.params['doubleSpinBox_2']['value'])
        
        acceptance = [float(i) for i in self.params['lineEdit_9']['value'].split(',')]
        if len(acceptance) == 1 and len(codes) > 1:
            acceptance = acceptance * len(codes)
        acceptance = np.array(acceptance)

        #kernels variables
        function = self.params['comboBox']['value']
        supports = [float(i) for i in self.params['lineEdit_5']['value'].split(',')]
        azms = [float(i) for i in self.params['lineEdit_2']['value'].split(',')]
        dips = [float(i) for i in self.params['lineEdit_3']['value'].split(',')]
        rakes = [float(i) for i in self.params['lineEdit_4']['value'].split(',')]
        nuggets = [float(i) for i in self.params['lineEdit_6']['value'].split(',')]
        major_med = [float(i) for i in self.params['lineEdit_7']['value'].split(',')]
        major_min = [float(i) for i in self.params['lineEdit_8']['value'].split(',')]
        
        kernels_par = {
        'supports':supports, 
        'azms':azms, 
        'dips':dips, 
        'rakes':rakes, 
        'nuggets':nuggets, 
        'major_meds':major_med, 
        'major_mins':major_min}

        for key in kernels_par:
            if len(kernels_par[key]) == 1 and len(codes) > 1:
                kernels_par[key] = kernels_par[key] * len(codes)

        kernels_par['function'] = [function] * len(codes)

        print(kernels_par)

        #getting coordinates
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

        #Interpolating variables
        a2g_grid = helpers.ar2gemsgrid_to_ar2gasgrid(tg_grid_name, tg_region_name)
        print('Computing results by RBF using a {} kernel...'.format(function))
        interpolated_variables = {}
        
        for idx, v in enumerate(variables):
            rt = int(codes[idx])
            interpolated_variables[rt] = {}
            print('Interpolating RT {}'.format(rt))

            if kernels_par['supports'][idx] == 0:
                mask = v < 0
                kernels_par['supports'][idx] = helpers.max_dist(x, y, z, mask, a2g_grid)
                print('Estimated support is {}'.format(kernels_par['supports'][idx]))

            rbf = RBF(x, y, z, v, function=kernels_par['function'][idx], nugget=kernels_par['nuggets'][idx], support=kernels_par['supports'][idx], major_med=kernels_par['major_meds'][idx], major_min=kernels_par['major_mins'][idx], azimuth=kernels_par['azms'][idx], dip=kernels_par['dips'][idx], rake=kernels_par['rakes'][idx])
            
            dc_param = dc_parametrization(v, dcf_param, f_min_o, f_min_p)
            for t in dc_param:
                print('Working on {}...'.format(t))
                
                rbf.train(dc_param[t])
                results = rbf.predict(a2g_grid)
                interpolated_variables[rt][t] = results

                if keep_variables == '1':
                    prop_name = 'interpolated_'+tg_prop_name+'_'+dcf_param+'_'+t+'_'+str(rt)
                    sgems.set_property(tg_grid_name, prop_name, results.tolist())

        print('Done!')

        #Generating geological models
        print('Building geomodel...')
        #getting closest node to each sample
        ids = [sgems.get_closest_nodeid(tg_grid_name, xi, yi, zi) for xi, yi, zi in zip(x,y,z)]
        rt_prop = sd_to_cat(variables, codes)
        build_geomodels(interpolated_variables, codes, tg_grid_name, tg_prop_name, ids, acceptance, rt_prop)

        print('Done!')

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "kernel_factor_uncertainty" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["kernel_factor_uncertainty"] #aqui vai o nome do plugin
