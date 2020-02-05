#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np

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

def entropy(prob_list):

	if True in np.isnan(prob_list):
		return float("nan")

	else:
		return -sum([prob*np.log2(prob) for prob in prob_list])

#################################################################################################

class block_transform: #aqui vai o nome do plugin
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
        grid_name = self.params['gridselectorbasic']['value']
        prop_names = self.params['orderedpropertyselector']['value'].split(';')
        gamma = float(self.params['doubleSpinBox']['value'])
        
        props_values = []
        for v in prop_names:
            props_values.append(sgems.get_property(grid_name, v))
        props_values = np.array(props_values)
        
        #calculating probs
        print('Calculating probabilities...')
        probs_matrix = np.array([sofmax_transformation(sds, gamma) for sds in props_values.T])
        probs_matrix = probs_matrix.T
        
        for i, p in enumerate(probs_matrix):
           sgems.set_property(grid_name, prop_names[i]+'_gamma_'+str(gamma), p.tolist())
           
        #calculating entropy
        print('Calculating entropy...')
        entropies = [entropy(probs) for probs in probs_matrix.T]
        entropies = (entropies - min(entropies))/(max(entropies) - min(entropies))
        sgems.set_property(grid_name, 'entropy_gamma_'+str(gamma), entropies.tolist())
        
        print('Done!')
        
        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "block_transform" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["block_transform"] #aqui vai o nome do plugin