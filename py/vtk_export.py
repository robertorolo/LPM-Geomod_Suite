#!/bin/python
#################################################################################################

#importe aqui os pacotes necessarios
import sgems
import numpy as np
from pyevtk.hl import pointsToVTK, gridToVTK

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

def grid_vertices(grid_name):
    info = sgems.get_grid_info(grid_name)
    nx, ny, nz = int(info['num_cells'][0]), int(info['num_cells'][1]), int(info['num_cells'][2]) 
    sx, sy, sz = float(info['dimension'][0]), float(info['dimension'][1]), float(info['dimension'][2]) 
    ox, oy, oz = float(info['origin'][0]), float(info['origin'][1]), float(info['origin'][2])
    X = np.array([(ox - sx/2 + x * sx) for x in range(nx+1)])
    Y = np.array([(oy - sy/2 + x * sy) for x in range(ny+1)])
    Z = np.array([(oz - sz/2 + x * sz) for x in range(nz+1)])

    return X, Y, Z

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
        pts_grid = self.params['gridselectorbasic']['value']
        pts_props = self.params['orderedpropertyselector']['value'].split(';')
        grids_grid = self.params['gridselectorbasic_2']['value']
        grids_props = self.params['orderedpropertyselector_2']['value'].split(';')
        flname_pts = self.params['filechooser']['value']
        flname_grid = self.params['filechooser_2']['value']

        if len(pts_props) != 0:
            x, y, z = np.array(sgems.get_X(pts_grid)), np.array(sgems.get_Y(pts_grid)), np.array(sgems.get_Z(pts_grid))
            vardict = {}
            for v in pts_props:
                values = np.array(sgems.get_property(pts_grid, v))
                vardict[v] = values
            pointsToVTK(flname_pts, x, y, z, data=vardict)
            print('Point set saved to {}'.format(flname_pts))

        if len(grids_props) !=0:
            X, Y, Z = grid_vertices(grids_grid)
            vardict = {}
            for v in grids_props:
                values = np.array(sgems.get_property(grids_grid, v))
                vardict[v] = values
            gridToVTK(flname_grid, X, Y, Z, cellData=vardict)
            print('Grid saved to {}'.format(flname_grid))

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "vtk_export" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["vtk_export"] #aqui vai o nome do plugin