#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
import math

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

def autogrid(x, y, z, sx, sy, sz, buffer=0):
	if z is None:
		nz = 1
		oz = 0
		max_z = 0
	else:
		oz = min(z) - buffer #+ sz/2
		max_z = max(z) + buffer + sz/2
		nz = math.ceil((max_z - oz)/sz)

	ox = min(x) - buffer #+ sx/2
	oy = min(y) - buffer #+ sy/2
	max_x = max(x) + buffer + sx/2
	max_y = max(y) + buffer + sy/2
	nx = math.ceil((max_x - ox)/(sx))
	ny = math.ceil((max_y - oy)/(sy))
	

	return {'ox':ox,'oy':oy,'oz':oz,'sx':sx,'sy':sy,'sz':sz,'nx':nx,'ny':ny,'nz':nz}

#################################################################################################

class auto_grid: #aqui vai o nome do plugin
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
        point_set_name = self.params['gridselectorbasic']['value']
        new_grid_name = self.params['lineEdit']['value']
        buffer = float(self.params['doubleSpinBox_4']['value'])
        sx, sy, sz = float(self.params['doubleSpinBox']['value']), float(self.params['doubleSpinBox_2']['value']), float(self.params['doubleSpinBox_3']['value'])

        x, y, z = np.array(sgems.get_X(point_set_name)), np.array(sgems.get_Y(point_set_name)), np.array(sgems.get_Z(point_set_name))
        grid_dic = autogrid(x, y, z, sx, sy, sz, buffer)

        sgems.execute('NewCartesianGrid  {}::{}::{}::{}::{}::{}::{}::{}::{}::{}::0,00'.format(new_grid_name, grid_dic['nx'], grid_dic['ny'], grid_dic['nz'], grid_dic['sx'], grid_dic['sy'], grid_dic['sz'], grid_dic['ox'], grid_dic['oy'], grid_dic['oz']))

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "auto_grid" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["auto_grid"] #aqui vai o nome do plugin