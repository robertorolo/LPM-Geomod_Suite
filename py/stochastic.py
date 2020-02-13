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
        
        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "stochastic" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["stochastic"] #aqui vai o nome do plugin