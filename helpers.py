import ar2gas
import sgems
from itertools import product
import numpy as np
from scipy.interpolate import NearestNDInterpolator

#################################################for an eventual use#################################################

def anis_dist(p0, p1, r1, r2, alpha, beta, gamma):
    

    if alpha >= 0 and alpha <=270:
        alpha = np.deg2rad(90-alpha)
    else:
        alpha = np.deg2rad(450-alpha)
    
    beta, gamma = np.deg2rad(beta), np.deg2rad(gamma)
    
    R = np.array([
    [np.cos(beta)*np.cos(alpha), np.cos(beta)*np.sin(alpha), -np.sin(beta)],
    [1/r1 * (-np.cos(gamma)*np.sin(alpha)+np.sin(gamma)*np.sin(beta)*np.cos(alpha)), 1/r1 * (np.cos(gamma)*np.cos(alpha)+np.sin(gamma)*np.sin(beta)*np.sin(alpha)), 1/r1 * (np.sin(gamma)*np.cos(beta))],
    [1/r2 * (np.sin(gamma)*np.sin(alpha)+np.cos(gamma)*np.cos(beta)*np.cos(alpha)), 1/r2 * (-np.sin(gamma)*np.cos(alpha)+np.cos(gamma)*np.sin(beta)*np.sin(alpha)), 1/r2 * (np.cos(gamma)*np.cos(beta))]
    ])
    
    p1_minus_p0 = np.array([
        [p1[0] - p0[0]],
        [p1[1] - p0[1]],
        [p1[2] - p0[2]]
    ])
        
    o1 = np.dot(p1_minus_p0.T, R.T)
    o2 = np.dot(o1,R)
    o3 = np.dot(o2, p1_minus_p0)
    
    return np.sqrt(o3)[0][0]

def spherical(h, a):
    if h < a:
        variance = 1-(3/2)*(h/a)+(1/2)*(h/a)**3
    else:
        variance = 1
    return variance

def gaussian(h, a):
    variance = np.exp(-(h/a)**2)
    return variance

def exponential(h, a):
    variance = np.exp(-h/a)
    return variance

#################################################for an eventual use#################################################


def modelfile_to_ar2gasmodel(path):
    f = open(path, "r")
    lines = f.readlines()
    cov_dict = {}
    for idx, l in enumerate(lines):
        if len(l) == 2:
            rt = int(l)
            model_structs = []
            nugget = float(lines[idx+1].split()[1].split('"')[1])
            model_structs.append(ar2gas.compute.Covariance.nugget(nugget))
            n_struct = int(lines[idx+1].split()[2].split('"')[1])
            for struc in range(n_struct):
                line = idx+struc*4+2
                contribution = float(lines[idx+2].split()[1].split('"')[1])
                struct_type = lines[idx+2].split()[2].split('"')[1]
                r1 = float(lines[line+1].split()[1].split('"')[1])
                r2 = float(lines[line+1].split()[2].split('"')[1])
                r3 = float(lines[line+1].split()[3].split('"')[1])
                a1 = float(lines[line+2].split()[1].split('"')[1])
                a2 = float(lines[line+2].split()[1].split('"')[1])
                a3 = float(lines[line+2].split()[1].split('"')[1])
                if struct_type=='Spherical':
                    struc_model = ar2gas.compute.Covariance.spherical(contribution, r1, r2, r3, a1, a2, a3)
                    model_structs.append(struc_model)
                if struct_type=='Exponential':
                    struc_model = ar2gas.compute.Covariance.exponential(contribution, r1, r2, r3, a1, a2, a3)
                    model_structs.append(struc_model)
                if struct_type=='Gaussian':
                    struc_model = ar2gas.compute.Covariance.gaussian(contribution, r1, r2, r3, a1, a2, a3)
                    model_structs.append(struc_model)
            cov_dict[rt] = model_structs
    return cov_dict

def ar2gemsvarwidget_to_ar2gascovariance(p):
    cov_lst = []
    n_var = int(p['indicator_regionalization_input']['number_of_indicator_group'])

    for i in range(n_var):
        indicator_group = "Indicator_group_" + str(i + 1)
        n_struct = int(p['indicator_regionalization_input'][indicator_group]['Covariance_input']['structures_count'])
        cov = []
        
        for j in range(n_struct):
            Structure = "Structure_" + str(j + 1)
            cov_type = p['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['type']

            cont = p['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['contribution']

            if cov_type == 'Nugget Covariance':
                nugget = ar2gas.compute.Covariance.nugget(float(cont))
                cov.append(nugget)

            else:
                range1 = p['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['ranges']['range1']
                range2 = p['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['ranges']['range2']
                range3 = p['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['ranges']['range3']

                rake = p['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['angles']['rake']
                dip = p['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['angles']['dip']
                azimuth = p['indicator_regionalization_input'][indicator_group]['Covariance_input'][Structure]['Two_point_model']['angles']['azimuth']

                if cov_type == 'Spherical Covariance':
                    struct = ar2gas.compute.Covariance.spherical(float(cont), float(range1), float(range2), float(range3), float(azimuth), float(dip), float(rake))
                    cov.append(struct)
                if cov_type == 'Exponential Covariance':
                    struct = ar2gas.compute.Covariance.exponential(float(cont), float(range1), float(range2), float(range3), float(azimuth), float(dip), float(rake))
                    cov.append(struct)
                if cov_type == 'Gaussian Covariance':
                    struct = ar2gas.compute.Covariance.gaussian(float(cont), float(range1), float(range2), float(range3), float(azimuth), float(dip), float(rake))
                    cov.append(struct)

        cov_lst.append(cov)

    return cov_lst

def ar2gemsgrid_to_ar2gasgrid(grid_name, region_name):
    info = sgems.get_grid_info(grid_name)
    grid = ar2gas.data.CartesianGrid(int(info['num_cells'][0]), int(info['num_cells'][1]), int(info['num_cells'][2]), 
                                         info['dimension'][0], info['dimension'][1], info['dimension'][2], 
                                         info['origin'][0], info['origin'][1], info['origin'][2])
    if region_name != '':
        region = np.array(sgems.get_region(grid_name, region_name))
        mask = region == 1
        #active_indexes = np.array([i for i in range(len(mask))])[mask]
        #grid = grid.grid_mask(mask.tolist())
        grid = ar2gas.data.MaskedGrid(grid, mask.tolist())
        #grid = ar2gas.data.MaskedGrid(grid, active_indexes.tolist())
        
    return grid

def ar2gasgrid_to_ar2gems(grid_name, grid):
    sgems.execute('NewCartesianGrid  {}::{}::{}::{}::{}::{}::{}::{}::{}::{}::0,00'.format(grid_name,
                                                                                          grid.dim()[0], grid.dim()[1], grid.dim()[2],
                                                                                          grid.cell_size()[0], grid.cell_size()[1], grid.cell_size()[2],
                                                                                          grid.origin()[0], grid.origin()[1], grid.origin()[2]))


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

def downscale_properties(grid, props, fx, fy, fz):
    if hasattr(grid, 'mask'):
        original_size = grid.size_of_mask()
        downscaled_grid = ar2gas.data.downscale_masked_grid(grid, fx, fy, fz)
        downscaled_size = downscaled_grid.size_of_mask()
    else:
        original_size = grid.size()
        downscaled_grid = ar2gas.data.downscale_cartesian_grid(grid, fx, fy, fz)
        downscaled_size = downscaled_grid.size()
    downscaled_props = []
    grid = ar2gas.data.CartesianGrid(grid.dim()[0], grid.dim()[1], grid.dim()[2], 1, 1, 1, 0, 0, 0) #bug fixing
    for p in props:
        downscaled_prop = np.zeros(downscaled_size)
        for i in range(original_size):
            value = p[i]
            ijk = grid.ijk(i)
            di = [i for i in range(ijk[0]*fx, ijk[0]*fx+fx)]
            dj = [j for j in range(ijk[1]*fy, ijk[1]*fy+fy)]
            dk = [k for k in range(ijk[2]*fz, ijk[2]*fz+fz)]
            
            for i, j, k in product(di, dj, dk):
                n = ijk_in_n(downscaled_grid, i, j, k)
                downscaled_prop[n] = value
        downscaled_props.append(downscaled_prop)
            
    return downscaled_grid, downscaled_props

def nn(x, y, z, var, grid):
    nan_mask = np.isfinite(var)
    points_array = np.vstack((x,y,z)).T
    knn = NearestNDInterpolator(points_array[nan_mask], var[nan_mask])
    grids_points_array = grid.locations()
    results = knn(grids_points_array)

    return results