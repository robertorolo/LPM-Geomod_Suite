#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import helpers
import ar2gas
from sklearn.metrics import confusion_matrix

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
        varg_path = path+'_varg_reproduction'
        con_mat_path = path+'_confusion_matrix'
        rt_grid_name = self.params['propertyselectornoregion']['grid']
        rt_prop_name = self.params['propertyselectornoregion']['property']
        n_var = int(self.params['indicator_regionalization_input']['number_of_indicator_group'])

        #getting codes
        rt_prop = np.array(sgems.get_property(rt_grid_name, rt_prop_name))
        nan_filter = np.isfinite(rt_prop)
        rt_prop_filterd = rt_prop[nan_filter]
        codes = np.unique(rt_prop_filterd)
        
        #calculating reals proportions
        reals = []
        for v in props_names:
            values = np.array(sgems.get_property(grid_name, v))
            reals.append(values)

        reals_props = prop_reals(reals, codes)
        nan_mask = np.isnan(reals[0])

        #calculating proportions
        print('Calculating proportions...')
        values = np.array(sgems.get_property(rt_grid_name, rt_prop_name))
        nan_mask_pts = np.isfinite(values)
        values = values[nan_mask_pts]
        x, y, z = np.array(sgems.get_X(rt_grid_name))[nan_mask_pts], np.array(sgems.get_Y(rt_grid_name))[nan_mask_pts], np.array(sgems.get_Z(rt_grid_name))[nan_mask_pts]
        grid = helpers.ar2gemsgrid_to_ar2gasgrid(grid_name, '')
        target = helpers.nn(x, y, z, values, grid)
        sgems.set_property(grid_name, 'target', target.tolist())
        #adding nan to target 
        for idx, i in enumerate(nan_mask):
            if i is True:
                target[idx] = float('nan')
                
        cat_dict = {}
        w = np.array([1/len(target)] * len(target))
        for c in codes:
            mask = target == c
            height = w[mask].sum()
            cat_dict[int(c)] = height
        
        #calculating variograms
        #experimental variograms
        print('Calculating variogrmas...')
        #getting variograms
        print('Getting variogram models')
        use_model_file = self.params['checkBox']['value'] 
        if use_model_file == '1':
            var_path = self.params['filechooser_2']['value']
            variograms = helpers.modelfile_to_ar2gasmodel(var_path)
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
        
        exp_vars_dict = {'ew':{}, 'ns':{}, 'z':{}}
        var_model_dict = {'ew':{}, 'ns':{}, 'z':{}}
        sx, sy, sz = grid.dim()[0], grid.dim()[1], grid.dim()[2]
        rangeinx, rangeiny, rangeinz = [i*sx for i in range(0,nlags+1)], [i*sy for i in range(0,nlags+1)], [i*sx for i in range(0,nlags+1)]

        print('Calculating experimental variograms...')
        for c in codes:
 
            indicators = np.array(reals) == c
            indicators = np.where(indicators==True, 1, 0)

            exp_vars_ew = [ar2gas.compute.compute_variogram(grid, data, (1, 0, 0), nlags, (0, 0, 0), 0) for data in indicators]
            exp_vars_dict['ew'][c] = exp_vars_ew

            exp_vars_ns = [ar2gas.compute.compute_variogram(grid, data, (0, 1, 0), nlags, (0, 0, 0), 0) for data in indicators]
            exp_vars_dict['ns'][c] = exp_vars_ns

            exp_vars_z = [ar2gas.compute.compute_variogram(grid, data, (0, 0, 1), nlags, (0, 0, 0), 0) for data in indicators]
            exp_vars_dict['z'][c] = exp_vars_z

        print('Calculating variogram models...')
        if np.sum(np.array(list(variograms.keys()))) == 0:
            c =0 
            cov = ar2gas.compute.KrigingCovariance(1., variograms[int(c)])
            sill = cov.compute([0,0,0],[0,0,0])
            
            model_var_ew = [sill-cov.compute([0,0,0],[pt,0,0]) for pt in rangeinx]
            var_model_dict['ew'][c] = model_var_ew

            model_var_ns = [sill-cov.compute([0,0,0],[0,pt,0]) for pt in rangeiny]
            var_model_dict['ns'][c] = model_var_ns

            model_var_z = [sill-cov.compute([0,0,0],[0,0,pt]) for pt in rangeinz]
            var_model_dict['z'][c] = model_var_z
        
        else:
            for c in codes:
                cov = ar2gas.compute.KrigingCovariance(1., variograms[int(c)])
                sill = cov.compute([0,0,0],[0,0,0])
            
                model_var_ew = [sill-cov.compute([0,0,0],[pt,0,0]) for pt in rangeinx]
                var_model_dict['ew'][c] = model_var_ew

                model_var_ns = [sill-cov.compute([0,0,0],[0,pt,0]) for pt in rangeiny]
                var_model_dict['ns'][c] = model_var_ns

                model_var_z = [sill-cov.compute([0,0,0],[0,0,pt]) for pt in rangeinz]
                var_model_dict['z'][c] = model_var_z
                
        #back flag
        print('Getting closest node for all realizations...')
        ids = [sgems.get_closest_nodeid(grid_name, xi, yi, zi) for xi, yi, zi in zip(x,y,z)]
        reals_values = [np.array(r)[ids] for r in reals]
        cms = [confusion_matrix(values, pred) for pred in reals_values]
        sum_ew = np.sum(cms, axis=0)
        final_cm = sum_ew / sum_ew.astype(np.float).sum(axis=1)

        print('Validation script saved at seleted folder!')

        script = '''
import numpy as np
import matplotlib.pyplot as plt

cat_dict = {}
reals_props = {}
colors = ['r','g','b','k','y','m','c']

def cat_plot(cat_dict, reals_props):
    #plotting target declustered histogram
    n = len(list(cat_dict.keys()))
    plt.bar(cat_dict.keys(), cat_dict.values(), color=colors[:n])
    plt.ylabel('Proportion')
    plt.xlabel('Categories')
    plt.xticks(list(cat_dict.keys()))
    #plotting realizations boxplots
    plt.boxplot(reals_props.values(), positions=list(cat_dict.keys()))
    plt.savefig('{}')

cat_plot(cat_dict, reals_props)

var_exp = {}
var_model = {}
model_keys = np.array(list(cat_dict.keys()))
ranges = [{}, {}, {}]
codes = list(cat_dict.keys())

def plt_vargs(var_exp, var_model):
    for c in codes:
        flname = 'C:/Users/robertomr/Documents/Testes/varg_reproduction'+'_'+str(c)
        fig, axes = plt.subplots(1,3, constrained_layout=True, figsize=(15,5))
        for idx, d in enumerate(['ew', 'ns', 'z']):
            #if np.all(model_keys):	
			#	axes[idx].plot(ranges[idx], var_model[d][0], color='red')
			#else:
			#	axes[idx].plot(ranges[idx], var_model[d][c], color='red')
            for r in var_exp[d][c]:
                axes[idx].plot(ranges[idx], r, color='gray')
                axes[idx].set_xlabel('Lag distance (m)')
                axes[idx].set_ylabel('Variance')
                axes[idx].set_title('Direction '+str(d))
                axes[idx].grid(True)
            
            if np.sum(model_keys) == 0:	
                axes[idx].plot(ranges[idx], var_model[d][0], color='red')
            else:
                axes[idx].plot(ranges[idx], var_model[d][c], color='red')
		#fig.title('Variogram '+str(int(c)))
        fig.savefig(flname)
		
plt_vargs(var_exp, var_model)

cm = {}

import seaborn as sns

plt.figure()
sns_plot = sns.heatmap(cm, annot=True, vmin=0.0, vmax=1.0, fmt='.2f')
plt.yticks(np.arange(len(codes))+0.5, labels=codes)
plt.xticks(np.arange(len(codes))+0.5, labels=codes)
plt.xlabel('Predicted')
plt.ylabel('Actual')
figure = sns_plot.get_figure()
figure.savefig('{}')

        '''.format(cat_dict, reals_props, hist_path, exp_vars_dict, var_model_dict, rangeinx, rangeiny, rangeinz, np.array2string(final_cm, separator=', '), con_mat_path)

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
