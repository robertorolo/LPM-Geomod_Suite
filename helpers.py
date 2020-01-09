import ar2gas

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
