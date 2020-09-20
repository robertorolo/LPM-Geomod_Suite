#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import helpers
import ar2gas
import time
from scipy.stats import norm
from sklearn.preprocessing import MinMaxScaler

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

def sample_negative(dists, pool_list, r):
    neg_f = dists <= 0
    neg_dists = dists[neg_f]
    neg_cats = pool_list[neg_f]
    neg_sorted_indexes = np.argsort(neg_dists)
    neg_dists_sorted = np.sort(neg_dists)
    how_many = len(neg_dists_sorted)

    if len(neg_dists_sorted) == 0:
        pos_f = dists >= 0
        pos_dists = dists[pos_f]
        pos_cats = pool_list[pos_f]
        cat = pos_cats[np.argmin(pos_dists)]

    # if there is only one negative distance
    if how_many == 1:
        if r < 0:
            cat = neg_cats[neg_sorted_indexes[0]]
    # if there is more than one negative distance
    else:
        for idx, i in enumerate(neg_dists_sorted):
            if idx == 0:
                if (r < neg_dists_sorted[idx+1]):
                    cat = neg_cats[neg_sorted_indexes[idx]]
            if idx == how_many-1:
                if (r > i) & (r < 0):
                    cat = neg_cats[neg_sorted_indexes[idx]]
            else:
                if (r > i) & (r < neg_dists_sorted[idx+1]):
                    cat = neg_cats[neg_sorted_indexes[idx]]
    
    return cat

def sample_positive(dists, pool_list, r):

    pos_f = dists >= 0
    pos_dists = dists[pos_f]
    pos_cats = pool_list[pos_f]
    pos_sorted_indexes = np.argsort(pos_dists)
    pos_dists_sorted = np.sort(pos_dists)
    how_many = len(pos_dists_sorted)

    if len(pos_dists_sorted) == 0:
        neg_f = dists <= 0
        neg_dists = dists[neg_f]
        neg_cats = pool_list[neg_f]
        cat = neg_cats[np.argmax(neg_dists)]

    # if there is only one positive distance
    if how_many == 1:
        if r > 0:
            cat = pos_cats[pos_sorted_indexes[0]]
    # if there is more than one positive distance
    else:
        for idx, i in enumerate(pos_dists_sorted):
            if idx == 0:
                if (r > 0) & (r < i):
                    cat = pos_cats[pos_sorted_indexes[idx]]
            if idx == how_many-1:
                if (r > pos_dists_sorted[idx-1]):
                    cat = pos_cats[pos_sorted_indexes[idx]]
            else:
                if (r > i) & (r < pos_dists_sorted[idx+1]):
                    cat = pos_cats[pos_sorted_indexes[idx]]

    return cat

def cat_sampler(dists, pool_list, r):
    if r > 0:
        cat = sample_positive(dists, pool_list, r)
    else:
        cat = sample_negative(dists, pool_list, r)

    return cat

def new_approach(dists, pool_list, r):
    cat = float('nan')
    sorted_indexes = np.argsort(dists)
    dists_sorted = np.sort(dists)

    if len(dists_sorted) == 2:
        if r < dists_sorted[1]:
            cat = pool_list[sorted_indexes[1]]
        else:
            cat = pool_list[sorted_indexes[0]]

    return cat

def num_cat(interpolated_variables, tg_grid_name, tg_prop_name, bw):

    num_cats = []
    for i in interpolated_variables.T:
        cats_inside = (i > -bw) & (i < bw)
        if np.any(cats_inside):
            num = np.sum(cats_inside)
            num_cats.append(num)
        else:
            num_cats.append(0)

    sgems.set_property(tg_grid_name, tg_prop_name+'_num_cats_{}'.format(bw), num_cats)

def build_binary_geomodel(interpolated_variable, tg_grid_name, tg_prop_name, real, bw):
    geomodel = np.ones(len(interpolated_variable)) * float('nan')
    indexes = np.arange(len(geomodel))
    outside = interpolated_variable > bw
    inside = interpolated_variable <= bw
    u_zone = (interpolated_variable > -bw) & (interpolated_variable < bw)
    indexes_u_zone = indexes[u_zone]
    geomodel[outside] = 0
    geomodel[inside] = 1 
    
    for i in indexes_u_zone:
        if real[i] > interpolated_variable[i]:
            geomodel[i] = 1
        else:
            geomodel[i] = 0

    sgems.set_property(tg_grid_name, tg_prop_name, geomodel.tolist())

def build_geomodel(interpolated_variables, codes, tg_grid_name, tg_prop_name, real, bw):
    
    geomodel = []
    for idx, i in enumerate(interpolated_variables.T):
        cats_inside = (i > -bw) & (i < bw)
        if np.sum(cats_inside) > 1:
            pool_list = np.array(codes)[cats_inside]
            dists = i[cats_inside]
            #r_val = np.random.choice(pool_list)
            #r_val = cat_sampler(dists, pool_list, real[idx])
            r_val = new_approach(dists, pool_list, real[idx])
            geomodel.append(r_val)
        else:
            index = i.argmin(axis=0)
            geomodel.append(float(codes[index]))

    sgems.set_property(tg_grid_name, tg_prop_name, geomodel)

def build_perturb_geomodel(interpolated_variables, codes, tg_grid_name, tg_prop_name, real, bw):
    
    geomodel = []
    for idx, i in enumerate(interpolated_variables.T):
        dists_inside = (i > -bw) & (i < bw)
        if np.sum(dists_inside) > 1:
            i[dists_inside] =  i[dists_inside]
            index = i.argmin(axis=0)
            #geomodel.append(float(codes[index]))
            geomodel.append(float('nan'))
        else:
            index = i.argmin(axis=0)
            geomodel.append(float(codes[index]))

    sgems.set_property(tg_grid_name, tg_prop_name, geomodel)

#################################################################################################

class boundary_simulation: #aqui vai o nome do plugin
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
        sd_grid_name = self.params['gridselectorbasic']['value']
        sd_var_names = self.params['orderedpropertyselector']['value']

        sim_grid_name = self.params['gridselectorbasic_2']['value']
        sim_var_names = self.params['orderedpropertyselector_2']['value']
        
        tg_grid_name = sim_grid_name
        tg_prop_name = self.params['lineEdit']['value'] 

        sd_var_names_lst = sd_var_names.split(';')
        sim_var_names_lst = sim_var_names.split(';')

        bw = float(self.params['doubleSpinBox']['value'])

        codes = [int(i.split('_')[-1]) for i in sd_var_names_lst]

        sd_matrix = np.array([sgems.get_property(sd_grid_name, p) for p in sd_var_names_lst])
        reals = np.array([sgems.get_property(sim_grid_name, p) for p in sim_var_names_lst])
        norm_reals = np.array([norm.cdf(lst) for lst in reals])
        scaler = MinMaxScaler(feature_range=(-bw, bw))
        #scaler = MinMaxScaler(feature_range=(-0.5, 0.5))
        reals_scaled = scaler.fit_transform(norm_reals)

        num_cat(sd_matrix, tg_grid_name, tg_prop_name, bw)
        
        for i, r in enumerate(reals_scaled):
            
            n_tg_prop_name = tg_prop_name+'_'+str(i)
            #build_geomodel(sd_matrix, codes, tg_grid_name, n_tg_prop_name, r, bw)
            build_perturb_geomodel(sd_matrix, codes, tg_grid_name, n_tg_prop_name, r, bw)
            
            #for l, var in enumerate(sd_matrix):
            #    n_tg_prop_name = tg_prop_name+'_'+str(codes[l])+'_'+str(i)
            #    build_binary_geomodel(var, tg_grid_name, n_tg_prop_name, r, bw)

        print('Finished!')

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "boundary_simulation" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["boundary_simulation"] #aqui vai o nome do plugin
