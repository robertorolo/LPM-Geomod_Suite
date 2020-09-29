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

def cat_random_sample(prob_list, u):

    position = 0
    probs = []
    for idx, prob in enumerate(prob_list):
        probs.append(prob)
        acc = sum(probs)
        if u <= acc:
            position = idx
            break

    return position

def prob_cat_sampler(dists, gamma, pool_list, r):
    prob_lst = sofmax_transformation(dists, gamma)
    a = cat_random_sample(prob_lst, r)
    cat = pool_list[a]

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

def build_geomodel(interpolated_variables, gamma, codes, tg_grid_name, tg_prop_name, real, bw):
    
    geomodel = []
    for idx, i in enumerate(interpolated_variables.T):
        cats_inside = (i > -bw) & (i < bw)
        if np.sum(cats_inside) > 1:
            pool_list = np.array(codes)[cats_inside]
            dists = i[cats_inside]
            #r_val = np.random.choice(pool_list)
            r_val = prob_cat_sampler(dists, gamma, pool_list, real[idx])
            geomodel.append(r_val)
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

        gamma = float(self.params['doubleSpinBox_2']['value'])

        sd_matrix = np.array([sgems.get_property(sd_grid_name, p) for p in sd_var_names_lst])
        if gamma == 0:
            gamma = np.max(sd_matrix)
        reals = np.array([sgems.get_property(sim_grid_name, p) for p in sim_var_names_lst])
        norm_reals = np.array([norm.cdf(lst) for lst in reals])
        #scaler = MinMaxScaler(feature_range=(-bw, bw))
        #reals_scaled = scaler.fit_transform(norm_reals)
        reals_scaled = norm_reals

        num_cat(sd_matrix, tg_grid_name, tg_prop_name, bw)
        
        for i, r in enumerate(reals_scaled):
            print('Working on real {}...'.format(i+1))

            if len(sd_matrix) > 1:
            
                n_tg_prop_name = tg_prop_name+'_'+str(i)
                build_geomodel(sd_matrix, gamma, codes, tg_grid_name, n_tg_prop_name, r, bw)
            
            else:
            
                for l, var in enumerate(sd_matrix):
                    n_tg_prop_name = tg_prop_name+'_'+str(codes[l])+'_'+str(i)
                    build_binary_geomodel(var, tg_grid_name, n_tg_prop_name, r, bw)

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
