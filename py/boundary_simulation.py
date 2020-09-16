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

def closest_cat(dists, pool_list, r):
    dists = np.absolute(np.array(dists) - r)
    print(r, dists, pool_list)
    index = dists.argmin(axis=0)
    cat = pool_list[index]
    print(cat)
    
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


def build_geomodel(interpolated_variables, codes, tg_grid_name, tg_prop_name, real, bw):
    
    geomodel = []
    for idx, i in enumerate(interpolated_variables.T):
        cats_inside = (i > -bw) & (i < bw)
        if np.any(cats_inside):
            pool_list = np.array(codes)[cats_inside]
            dists = i[cats_inside]
            r_val = sample_cat(dists, pool_list, real[idx])
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

        sd_matrix = np.array([sgems.get_property(sd_grid_name, p) for p in sd_var_names_lst])
        reals = np.array([sgems.get_property(sim_grid_name, p) for p in sim_var_names_lst])
        norm_reals = np.array([norm.cdf(lst) for lst in reals])
        scaler = MinMaxScaler(feature_range=(-bw, bw))
        reals_scaled = scaler.fit_transform(norm_reals)

        num_cat(sd_matrix, tg_grid_name, tg_prop_name, bw)
        
        for i, r in enumerate(reals_scaled):
            
            n_tg_prop_name = tg_prop_name+'_'+str(i)
            build_geomodel(sd_matrix, codes, tg_grid_name, n_tg_prop_name, r, bw)

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
