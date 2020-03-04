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
  
def std_array(a, var):
    a = np.array(a)
    std_a = (a - np.nanmin(a))/(np.nanmax(a) - np.nanmin(a))
    std_a = std_a * var
    return std_a.tolist()

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
      
        print('Getting variables...')
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
        std_var = self.params['checkBox_2']['value']

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
        #adding nan to target 
        target[nan_mask]=float('nan')
        sgems.set_property(grid_name, 'target', target.tolist())
                
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
        use_model_file = self.params['checkBox_3']['value'] 
        if use_model_file == '1':
            var_path = self.params['filechooser_3']['value']
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
        
        dimz = grid.dim()[2]
        pts_exp_vars = 0
        if dimz > 1:
            dimension = 3
        else:
            dimension = 2
        plot_exp = self.params['checkBox']['value']
        if plot_exp == '1':
            exp_path = self.params['filechooser_2']['value'] 
            pts_exp_vars = helpers.read_exp_vars(exp_path, dimension)
        
        exp_vars_dict = {'ew':{}, 'ew_mean':{}, 'ns':{}, 'ns_mean':{}, 'z':{}, 'z_mean':{}}
        var_model_dict = {'ew':{}, 'ns':{}, 'z':{}}
        sx, sy, sz = grid.cell_size()[0], grid.cell_size()[1], grid.cell_size()[2]
        rangeinx, rangeiny, rangeinz = [i*sx for i in range(0,nlags+1)], [i*sy for i in range(0,nlags+1)], [i*sz for i in range(0,nlags+1)]

        print('Calculating experimental variograms...')
        for c in codes:
 
            indicators = np.array(reals) == c
            indicators = np.where(indicators==True, 1, 0)

            exp_vars_ew = [ar2gas.compute.compute_variogram(grid, data, (1, 0, 0), nlags, (0, 0, 0), 0) for data in indicators]
            mean_ew = np.mean(np.array(exp_vars_ew), axis=0)
            exp_vars_dict['ew'][c] = exp_vars_ew
            exp_vars_dict['ew_mean'][c] = mean_ew.tolist()

            exp_vars_ns = [ar2gas.compute.compute_variogram(grid, data, (0, 1, 0), nlags, (0, 0, 0), 0) for data in indicators]
            mean_ns = np.mean(np.array(exp_vars_ns), axis=0)
            exp_vars_dict['ns'][c] = exp_vars_ns
            exp_vars_dict['ns_mean'][c] = mean_ns.tolist()

            if dimz > 1:
            
                exp_vars_z = [ar2gas.compute.compute_variogram(grid, data, (0, 0, 1), nlags, (0, 0, 0), 0) for data in indicators]
                mean_z = np.mean(np.array(exp_vars_z), axis=0)
                exp_vars_dict['z'][c] = exp_vars_z
                exp_vars_dict['z_mean'][c] = mean_z.tolist()

        print('Calculating variogram models...')
        for c in codes:
            if np.sum(np.array(list(variograms.keys()))) == 0:
                c1 = 0
            else:
                c1 = c 
            cov = ar2gas.compute.KrigingCovariance(1., variograms[int(c1)])
            sill = cov.compute([0,0,0],[0,0,0])
        
            model_var_ew = [sill-cov.compute([0,0,0],[pt,0,0]) for pt in rangeinx]
            if std_var == '1':
                var_model_dict['ew'][c1] = std_array(model_var_ew, np.nanmax(exp_vars_dict['ew_mean'][c]))
            else:
                var_model_dict['ew'][c1] = model_var_ew

            model_var_ns = [sill-cov.compute([0,0,0],[0,pt,0]) for pt in rangeiny]
            if std_var == '1':
                var_model_dict['ns'][c1] = std_array(model_var_ns, np.nanmax(exp_vars_dict['ns_mean'][c]))
            else:
                var_model_dict['ns'][c1] = model_var_ns

            if dimz > 1:
            
                model_var_z = [sill-cov.compute([0,0,0],[0,0,pt]) for pt in rangeinz]
                if std_var == '1':
                    var_model_dict['z'][c1] = std_array(model_var_z, np.nanmax(exp_vars_dict['z_mean'][c]))
                else:
                    var_model_dict['z'][c1] = model_var_z
                
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

nan = float('nan')

plot_exp = int({})
pts_exp_var = {}

dimz = {}
if dimz > 1:
    dirs = ['ew', 'ns', 'z']
    dp = 3
else:
    dirs = ['ew', 'ns']
    dp = 2

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
flname = '{}'

def plt_vargs(var_exp, var_model, flname):
    for c in codes:
        flname_c = flname+'_'+str(c)
        fig, axes = plt.subplots(1,dp, constrained_layout=True, figsize=(15,5))
        for idx, d in enumerate(dirs):
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
            
            axes[idx].plot(ranges[idx], var_exp[d+str('_mean')][c], color='black')
            if len(var_model[d]) == 1:	
                axes[idx].plot(ranges[idx], var_model[d][0], color='red')
            else:
                axes[idx].plot(ranges[idx], var_model[d][c], color='red')
            
            if plot_exp == 1:
                axes[idx].plot(pts_exp_var[c][d]['x'], pts_exp_var[c][d]['y'], linestyle='dashed', marker='x', color='blue')
                
		#fig.title('Variogram '+str(int(c)))
        fig.savefig(flname_c)
		
plt_vargs(var_exp, var_model, flname)

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

        '''.format(plot_exp, pts_exp_vars, dimz, cat_dict, reals_props, hist_path, exp_vars_dict, var_model_dict, rangeinx, rangeiny, rangeinz, varg_path, np.array2string(final_cm, separator=', '), con_mat_path)

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
