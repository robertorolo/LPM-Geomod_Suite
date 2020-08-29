#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import helpers
import ar2gas
import time
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

#################################################################################################

class pfields: #aqui vai o nome do plugin
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
        prob_grid_name = self.params['gridselectorbasic']['value']
        prob_var_names = self.params['orderedpropertyselector']['value']

        sim_grid_name = self.params['gridselectorbasic_2']['value']
        sim_var_names = self.params['orderedpropertyselector_2']['value']
        
        tg_grid_name = sim_grid_name
        tg_prop_name = self.params['lineEdit']['value'] 

        prob_var_names_lst = prob_var_names.split(';')
        sim_var_names_lst = sim_var_names.split(';')

        
        codes = [int(i.split('_')[-3]) for i in prob_var_names_lst]

        probs_matrix = np.array([sgems.get_property(prob_grid_name, p) for p in prob_var_names_lst])
        reals = np.array([sgems.get_property(sim_grid_name, p) for p in sim_var_names_lst])
        norm_reals = [norm.cdf(lst) for lst in reals]

        #sampling cats
        print('Sampling categories using the p-fields and building realizations...')
        t1 = time.time()
        
        probs_matrix = probs_matrix.T
        for real_idx, r in enumerate(norm_reals):
            realization = []
            for idx, b in enumerate(np.array(r).T):
                if np.isnan(probs_matrix[idx]).any():
                    realization.append(float('nan'))
                else:
                    position = cat_random_sample(probs_matrix[idx], b)
                    realization.append(int(codes[position]))
            sgems.set_property(tg_grid_name, tg_prop_name+'_real_'+str(real_idx), realization)

        t2 = time.time()
        print('Took {} seconds'.format(round((t2-t1), 2)))

        print('Finished!')

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "pfields" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["pfields"] #aqui vai o nome do plugin
