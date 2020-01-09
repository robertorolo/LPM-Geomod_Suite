import ar2gas
import sgems
from itertools import product

def modelfile_to_ar2gasmodel(path):
    f = open(path, "r")
    l = f.readlines()
    model_structs = []
    nugget = float(l[0].split()[1].split('"')[1])
    model_structs.append(ar2gas.compute.Covariance.nugget(nugget))
    n_struct = int(l[0].split()[2].split('"')[1])
    for struc in range(n_struct):
        line = struc*4+1
        contribution = float(l[line].split()[1].split('"')[1])
        struct_type = l[line].split()[2].split('"')[1]
        r1 = float(l[line+1].split()[1].split('"')[1])
        r2 = float(l[line+1].split()[2].split('"')[1])
        r3 = float(l[line+1].split()[3].split('"')[1])
        a1 = float(l[line+2].split()[1].split('"')[1])
        a2 = float(l[line+2].split()[1].split('"')[1])
        a3 = float(l[line+2].split()[1].split('"')[1])
        if struct_type=='Spherical':
            struc_model = ar2gas.compute.Covariance.spherical(contribution, r1, r2, r3, a1, a2, a3)
            model_structs.append(struc_model)
        if struct_type=='Exponential':
            struc_model = ar2gas.compute.Covariance.exponential(contribution, r1, r2, r3, a1, a2, a3)
            model_structs.append(struc_model)
        if struct_type=='Gaussian':
            struc_model = ar2gas.compute.Covariance.gaussian(contribution, r1, r2, r3, a1, a2, a3)
            model_structs.append(struc_model)
    return model_structs

def ar2gemsgrid_to_ar2gasgrid(grid_name, region_name=None):
    info = sgems.get_grid_info(grid_name)
    
    if region_name is None:
        grid = ar2gas.data.CartesianGrid(region, 
                                         info['num_cells'][0], info['num_cells'][1], info['num_cells'][2], 
                                         info['dimension'][0], info['dimension'][1], info['dimension'][2], 
                                         info['origin'][0], info['origin'][1], info['origin'][2])
    else:
        region = sgems.get_region(grid_name, region)
        grid = ar2gas.data.MaskedGrid(grid, region)
        
    return grid

def ar2gasprop_to_ar2gems(grid, grid_name, prop, prop_name):
    sgems.execute('NewCartesianGrid  {}::{}::{}::{}::{}::{}::{}::{}::{}::{}::0,00'.format(grid_name,
                                                                                          grid.dim()[0], grid.dim()[1], grid.dim()[2],
                                                                                          grid.cell_size()[0], grid.cell_size()[1], grid.cell_size()[2],
                                                                                          grid.origin()[0], grid.origin()[1], grid.origin()[2]))
    sgems.set_property(grid_name, prop_name, prop)

def ijk_in_n(grid, i, j, k):
    dims = grid.dim()
    n = k*dims[0]*dims[1]+j*dims[0]+i
    return n

def downscale_property(grid, prop, fx, fy, fz):
    downscaled_grid = ar2gas.data.downscale_cartesian_grid(grid, fx, fy, fz)
    downscaled_prop = np.zeros(downscaled_grid.size())
    for i in range(grid.size()):
        value = prop[i]
        ijk = grid.ijk(i)
        di = [i for i in range(ijk[0]*fx, ijk[0]*fx+fx)]
        dj = [j for j in range(ijk[1]*fy, ijk[1]*fy+fy)]
        dk = [k for k in range(ijk[2]*fz, ijk[2]*fz+fz)]
        
        for i, j, k in product(di, dj, dk):
            n = ijk_in_n(downscaled_grid, i, j, k)
            downscaled_prop[n] = value
            
    return downscaled_grid, downscaled_prop.tolist()
