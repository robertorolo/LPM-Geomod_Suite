#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
from scipy.spatial import distance

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

def min_dist(point, points):
	dist_matrix = distance.cdist([point], points, "euclidean")
	return np.amin(dist_matrix)

#################################################################################################

class signed_distances_transform: #aqui vai o nome do plugin
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
        x, y, z = np.array(sgems.get_X(grid_name)), np.array(sgems.get_Y(grid_name)), np.array(sgems.get_Z(grid_name))
        coords_matrix = np.vstack((x,y,z)).T

	    #calculating signed distances
        nan_filter = np.isfinite(prop_values)
        unique_rts = np.unique(prop_values[nan_filter])

        for rt in unique_rts:
            print('calculating signed distances for rock type {}'.format(int(rt)))
            filter_0 = prop_values != rt
            filter_1 = prop_values == rt
            points_0 = coords_matrix[filter_0]
            points_1 = coords_matrix[filter_1]

            sd_prop = []
            for idx, pt in enumerate(prop_values):
                if np.isnan(pt):
                    sd_prop.append(float('nan'))

                else:
                    point = coords_matrix[idx]
                    if pt == rt:
                        sd_prop.append(-min_dist(point, points_0))
                    else:
                       sd_prop.append(min_dist(point, points_1)) 
        
            sgems.set_property(grid_name, 'signed_distances_rt_{}'.format(int(rt)), sd_prop)
        
        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "signed_distances_transform" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["signed_distances_transform"] #aqui vai o nome do plugin