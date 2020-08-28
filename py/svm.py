#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import helpers
import ar2gas
import time

from sklearn.preprocessing import MinMaxScaler
from sklearn.multiclass import OneVsOneClassifier
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV

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

def str_to_float(sting):
    splitted = sting.split('**')
    float_val = float(splitted[0])**float(splitted[1])
    return float_val

#################################################################################################

class SVM: #aqui vai o nome do plugin
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
        tg_grid_name = self.params['gridselector']['value']
        tg_region_name = self.params['gridselector']['region']
        tg_prop_name = self.params['lineEdit']['value']
        prop_grid_name = self.params['propertyselectornoregion']['grid'] 
        prop_name = self.params['propertyselectornoregion']['property']

        cmin = str_to_float(self.params['lineEdit_3']['value'])
        cmax = str_to_float(self.params['lineEdit_2']['value'])
        gamma_min = str_to_float(self.params['lineEdit_5']['value'])
        gamma_max = str_to_float(self.params['lineEdit_4']['value'])
        n = int(self.params['spinBox']['value'])

        n_folds = int(self.params['spinBox_2']['value'])

        values = np.array(sgems.get_property(prop_grid_name, prop_name))
        nan_filter = np.isfinite(values)
        y_val = values[nan_filter]

        x, y, z = np.array(sgems.get_X(prop_grid_name))[nan_filter], np.array(sgems.get_Y(prop_grid_name))[nan_filter], np.array(sgems.get_Z(prop_grid_name))[nan_filter]
        X = np.array([x, y, z]).T 

        grid = helpers.ar2gemsgrid_to_ar2gasgrid(tg_grid_name, tg_region_name)
        X_new = grid.locations()

        #scaling
        scaler = MinMaxScaler(feature_range=(0, 1))
        X = scaler.fit_transform(X)
        X_new = scaler.fit_transform(X_new)

        #grid search parameters
        C = np.linspace(cmin, cmax, n)
        gamma = np.linspace(gamma_min, gamma_max, n)
        param_grid = [
        {'estimator__C': C, 'estimator__gamma': gamma, 'estimator__kernel':["rbf"]*n},
        ]

        #clf
        print("Tunning parameters...")
        clf = OneVsOneClassifier(SVC())
        model_tunning = GridSearchCV(clf, param_grid=param_grid, cv=n_folds)
        model_tunning.fit(X, y_val)

        print('accuracy ', model_tunning.best_score_)
        print('parameters', model_tunning.best_params_)

        clf = model_tunning.best_estimator_

        print("Predicting...")
        clf.fit(X, y_val)
        results = clf.predict(X_new)

        sgems.set_property(tg_grid_name, tg_prop_name, results.tolist())
        print('Done!')

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "SVM" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["SVM"] #aqui vai o nome do plugin
