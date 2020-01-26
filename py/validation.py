#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import helpers

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
    props[int(c)] = cat_props
  return props

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
        scipt_path = path+'.py'
        hist_path = path+'_hist_reproduction'
        rt_grid_name = self.params['propertyselectornoregion']['grid']
        rt_prop_name = self.params['propertyselectornoregion']['property']

        #getting codes
        rt_prop = np.array(sgems.get_property(rt_grid_name, rt_prop_name))
        nan_filter = np.isfinite(rt_prop)
        rt_prop_filterd = rt_prop[nan_filter]
        codes = np.unique(rt_prop_filterd)
        my_codes_string = ','.join(map(str, codes.astype(int))) 

        #calculating target proportions
        x, y, z = np.array(sgems.get_X(rt_grid_name)), np.array(sgems.get_Y(rt_grid_name)), np.array(sgems.get_Z(rt_grid_name))
        values = np.array(sgems.get_property(rt_grid_name, rt_prop_name))
        grid = helpers.ar2gemsgrid_to_ar2gasgrid(grid_name, '')
        target = helpers.nn(x, y, z, values, grid)
        cat_dict = {}
        w = np.array([1/len(target)] * len(target))
        for c in codes:
            mask = target == c
            height = w[mask].sum()
            cat_dict[int(c)] = height

        reals = []
        for v in props_names:
            values = np.array(sgems.get_property(grid_name, v))
            reals.append(values)

        reals_props = prop_reals(reals, codes)

        script = '''
import numpy as np
import matplotlib.pyplot as plt

cat_dict = {}
reals_props = {}
cats = [{}]

def cat_plot(cat_dict, reals_props, cats):
    #plotting target declustered histogram
    plt.bar(cat_dict.keys(), cat_dict.values(), color='gray')
    plt.ylabel('Proportion')
    plt.xlabel('Categories')
    plt.xticks(cats)
    #plotting realizations boxplots
    plt.boxplot(reals_props.values(), positions=cats)
    plt.savefig('{}')

cat_plot(cat_dict, reals_props, cats)
        '''.format(cat_dict, reals_props, my_codes_string, hist_path)

        #writing script
        f = open(scipt_path, 'w')
        f.write(script)
        f.close()

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "validation" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["validation"] #aqui vai o nome do plugin