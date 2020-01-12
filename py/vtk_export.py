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

class vtk_export: #aqui vai o nome do plugin
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

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "vtk_export" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["vtk_export"] #aqui vai o nome do plugin