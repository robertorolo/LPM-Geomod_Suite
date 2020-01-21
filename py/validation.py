#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import helpers
import matplotlib.pyplot as plt

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

def prop_reals(reals, cats):
  props = {}
  for c in cats:
    cat_props = []
    for r in reals:
      cat_filter = np.array(r) == c
      real_prop = len(np.array(r)[cat_filter])/len(r)
      cat_props.append(real_prop)
    props[c] = cat_props
  return props

def cat_plot(target, weights, reals, path):
    #plotting target declustered histogram
    weights = np.array(weights) if weights is not None else np.array([1/len(target)] * len(target))
    cats = np.unique(target)
    cat_dict = {}
    for c in cats:
        mask = target == c
        height = weights[mask].sum()
        cat_dict[c] = height
    plt.figure(figsize=(5,8))
    plt.bar(cat_dict.keys(), cat_dict.values())
    plt.ylabel('Proportion')
    plt.xlabel('Categories')
    plt.xticks(list(cat_dict.keys()))
    #plotting realizations boxplots
    reals_props = prop_reals(reals, cats)
    plt.boxplot(reals_props.values(), positions=cats)
    plt.savefig(path)

#################################################################################################

class validation: #aqui vai o nome do plugin
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
        props_names = self.params['orderedpropertyselector']['value'].split(';')
        nlags = int(self.params['spinBox']['value'])
        path = self.params['filechooser']['value']
        rt_grid_name = self.params['propertyselectornoregion']['grid']
        rt_prop_name = self.params['propertyselectornoregion']['property']

        #getting codes
        rt_prop = np.array(sgems.get_property(rt_grid_name, rt_prop_name))
        nan_filter = np.isfinite(rt_prop)
        rt_prop_filterd = rt_prop[nan_filter]
        codes = np.unique(rt_prop_filterd)

        #transforming grid to ar2gas grid
        a2g_grid = helpers.ar2gemsgrid_to_ar2gasgrid(grid_name, '')

        #getting declustered rt
        x, y, z = np.array(sgems.get_X(rt_grid_name))[nan_filter], np.array(sgems.get_Y(rt_grid_name))[nan_filter], np.array(sgems.get_Z(rt_grid_name))[nan_filter]
        nn_results = helpers.nn(x, y, z, rt_prop_filterd, a2g_grid)

        #getting reals props 
        realizations = []
        for v in props_names:
            values =  np.array(sgems.get_property(grid_name, v))
            realizations.append(values)

        #getting variograms
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
            if len(varg_lst) == 1:
                varg_lst = varg_lst * len(codes)
                variograms = {}
                variograms[0] = varg_lst[0]
            else:
                variograms = dict(zip(codes, varg_lst))

        #saving histogram
        cat_plot(nn_results, None, realizations, path)

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "validation" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["validation"] #aqui vai o nome do plugin