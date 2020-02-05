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

class block_transform: #aqui vai o nome do plugin
    def __init__(self):
        pass

#################################################################################################

    def initialize(self, params):
        self.params = params
        
        #imprimindo o dicionario de parametros
        print("dicionario de parametros: ", params)

        #executando a funcao exibe os valores do dicionario de parametros
        read_params(params) #para nao printar comente essa linha

        return True

#################################################################################################

    def execute(self):
      
        #aqui vai o codigo
        #getting variables

        return True

#################################################################################################

    def finalize(self):
        return True

    def name(self):
        return "block_transform" #aqui vai o nome do plugin

#################################################################################################

def get_plugins():
    return ["block_transform"] #aqui vai o nome do plugin