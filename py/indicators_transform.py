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

#################################################################################################

class indicators_transform: #aqui vai o nome do plugin
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
        grid_name = self.params['rt_prop']['grid']
        prop_name = self.params['rt_prop']['property']
        prop_values = sgems.get_property(grid_name, prop_name)
        prop_values = np.array(prop_values)

	#calculating indicators
        nan_filter = np.isfinite(prop_values)
        unique_rts = np.unique(prop_values[nan_filter])

        for rt in unique_rts:
            print('calculating indicators for rock type {}'.format(int(rt)))
            ind_prop = prop_values == rt
            ind_prop = np.where(ind_prop == True, 1., 0.)
            ind_prop[~nan_filter] = float('nan')
            indprop_name = 'indicators_rt_{}'.format(int(rt))
            sgems.set_property(grid_name, indprop_name, ind_prop.tolist())

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "indicators_transform" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["indicators_transform"] #aqui vai o nome do plugin