import numpy as np

def load(filename="pot.mtp"):
    f = open(filename)
    lines = f.read().splitlines()
    f.close()

    assert lines[0] == "MTP", "Can read only MTP format potentials"
    version            = lines[1].split("=")[1].strip()
    
    assert version == "1.1.0", "MTP file must have version \"1.1.0\""
    pot_desc           = lines[2].split("=")[1].strip()
    scaling            = float(lines[3].split("=")[1])
    species_count      = int(lines[4].split("=")[1])
    # potential_tag      = lines[5].split("=")[1].strip()
    rbasis_type        = lines[6].split("=")[1].strip()
    min_dist           = float(lines[7].split("=")[1])
    max_dist           = float(lines[8].split("=")[1])
    rb_size            = int(lines[9].split("=")[1])
    
    # associate rbasis_type  to the correct potnteital pointer
    # scaling pointer

    rb_vals = np.zeros(rb_size)
    


    radial_func_count  = int(lines[10].split("=")[1])

    # Radial coeffs initialization
    pairs_count = species_count ** 2; #number of species pairs
    n   = pairs_count * radial_func_count * rb_size # pairs_count*radial_func_count*(p_RadialBasis->rb_size)
    regression_coeffs = np.zeros(n)

    # inited = True;

    iline = 12
    for s1 in range(species_count):
        m = s1 * species_count
        for s2 in range(species_count):
            iline += 1
            n = (m + s2) * radial_func_count * rb_size
            for i in range(radial_func_count):
                p = n + (i * rb_size)
                list_t = lines[iline].strip("\t").strip("{").strip("}").split(",")
                iline += 1
                for j in range(rb_size):
                    regression_coeffs[p + j] = list_t[j]
    #

    alpha_moments_count = int(lines[iline].split("=")[1])
    alpha_index_basic_count = int(lines[iline + 1].split("=")[1])

    vtemp = lines[iline + 2].split("=")[1]
    vtemp = vtemp.strip(" {{").strip("}}").split("}, {")
    alpha_index_basic = np.zeros((alpha_index_basic_count, 4), int)
    radial_func_max = -1
    for i in range(alpha_index_basic_count):
        v = vtemp[i].split(",")
        for j in range(4):
            alpha_index_basic[i,j] = int(v[j])
        #
        alphaIB0 = alpha_index_basic[i, 0]
        if alphaIB0 > radial_func_max:
            radial_func_max = alphaIB0
        #
    #
    assert radial_func_max == (radial_func_count - 1), "Wrong number of radial functions specified"


    alpha_index_times_count = int(lines[iline + 3].split("=")[1])

    vtemp = lines[iline + 4].split("=")[1]
    vtemp = vtemp.strip(" {{").strip("}}").split("}, {")
    alpha_index_times = np.zeros((alpha_index_times_count, 4), int)
    for i in range(alpha_index_times_count):
        v = vtemp[i].split(",")
        for j in range(4):
            alpha_index_times[i,j] = int(v[j])
        #
    #


    alpha_scalar_moments = int(lines[iline + 5].split("=")[1])
    vtemp = lines[iline + 6].split("=")[1].strip().strip("{").strip("}").split(",")
    alpha_moment_mapping = [int(i) for i in vtemp]

    alpha_count = alpha_scalar_moments + 1;

    # species_coeffs
    vtemp1 = lines[iline + 7].split("=")[1].strip().strip("{").strip("}").split(",")
    # linear_coeffs = [int(i) for i in vtemp]
    # assert len(linear_coeffs) == species_count

    # moment_coeffs
    vtemp2 = lines[iline + 8].split("=")[1].strip().strip("{").strip("}").split(",")
    # linear_coeffs = [int(i) for i in vtemp]

    l1 = len(vtemp1)
    l  = l1 + len(vtemp2)
    linear_coeffs = np.zeros(l)
    for i in range(l1):
        linear_coeffs[i] = vtemp1[i]
    #
    j = 0
    for i in range(l1, l):
        linear_coeffs[i] = vtemp2[j]
        j += 1
    #
    #############
    # n = alpha_count - 1 + species_count
    # energy_cmpnts = np.zero(n)
    # forces_cmpnts = np.zeros(n * 3)
    # stress_cmpnts = (double(*)[3][3])malloc(n * sizeof(stress_cmpnts[0]));

    moment_vals = np.zeros(alpha_moments_count)
    # basis_vals = np.zeros(alpha_count)
    # site_energy_ders_wrt_moments_ = np.zeros(alpha_moments_count)
    #############
    radial_size = len(regression_coeffs)
    max_comp = species_count - 1

    regression_coeffs = [regression_coeffs[i] for i in range(len(regression_coeffs))]
    for e in linear_coeffs:
        regression_coeffs.append(e)
    #
    
    linear_mults = np.ones(alpha_scalar_moments)
    max_linear = 1e-10 * np.ones(alpha_scalar_moments)
#

