import numpy as np
from ase import neighborlist

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

    # rb_vals = np.zeros(rb_size)
    


    radial_func_count  = int(lines[10].split("=")[1])

    # Radial coeffs initialization
    pairs_count = species_count ** 2; #number of species pairs
    n   = pairs_count * radial_func_count * rb_size # pairs_count*radial_func_count*(p_RadialBasis->rb_size)
    # regression_coeffs = np.zeros(n)

    iline = 12
    # for s1 in range(species_count):
    #     m = s1 * species_count
    #     for s2 in range(species_count):
    #         iline += 1
    #         n = (m + s2) * radial_func_count * rb_size
    #         for i in range(radial_func_count):
    #             p = n + (i * rb_size)
    #             list_t = lines[iline].strip("\t").strip("{").strip("}").split(",")
    #             iline += 1
    #             for j in range(rb_size):
    #                 regression_coeffs[p + j] = list_t[j]
    # #

    regression_coeffs = np.zeros( (species_count, species_count, radial_func_count, rb_size) )
    for s1 in range(species_count):
        for s2 in range(species_count):
            iline += 1
            for i in range(radial_func_count):
                list_t = lines[iline].strip("\t").strip("{").strip("}").split(",")
                iline += 1
                for j in range(rb_size):
                    regression_coeffs[s1, s2, i, j] = list_t[j]
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

    # moment_vals = np.zeros(alpha_moments_count)
    # basis_vals = np.zeros(alpha_count)
    # site_energy_ders_wrt_moments_ = np.zeros(alpha_moments_count)
    #############
    # radial_size = len(regression_coeffs)
    # max_comp = species_count - 1

    # regression_coeffs = [regression_coeffs[i] for i in range(len(regression_coeffs))]
    # for e in linear_coeffs:
    #     regression_coeffs.append(e)
    # #
    
    # toReturn = [pot_desc, scaling, species_count, 
    #             rbasis_type, min_dist, max_dist, 
    #             rb_size, radial_func_count, regression_coeffs,
    #             alpha_moments_count, alpha_index_basic_count, 
    #             alpha_index_basic, alpha_index_times, 
    #             alpha_scalar_moments, linear_coeffs]
    #
    params = { "pot_desc":pot_desc,
            "scaling":scaling,
            "species_count":species_count,
            "rbasis_type":rbasis_type,
            "min_dist":min_dist,
            "max_dist":max_dist,
            "rb_size":rb_size,
            "radial_func_count":radial_func_count,
            "regression_coeffs":regression_coeffs,
            "alpha_moments_count":alpha_moments_count,
            "alpha_index_basic_count":alpha_index_basic_count,
            "alpha_index_basic":alpha_index_basic,
            "alpha_index_times":alpha_index_times,
            "alpha_scalar_moments":alpha_scalar_moments,
            "alpha_count":alpha_count,
            "linear_coeffs":linear_coeffs,
            "alpha_index_times_count":alpha_index_times_count,
            "alpha_moment_mapping":alpha_moment_mapping,
    }
    return params
#

def init_vecs(parameters):
    rb_vals = np.zeros(parameters['rb_size'])
    rb_ders = np.zeros(parameters["rb_size"])

    moment_vals = np.zeros(parameters['alpha_moments_count'])
    basis_vals = np.zeros(parameters['alpha_count']) 
    site_energy_ders_wrt_moments_ = np.zeros(parameters['alpha_moments_count'])
    max_dist = parameters['max_dist']
    min_dist = parameters['min_dist']
    mult = 2.0 / (max_dist - min_dist)
    
    max_alpha_index_basic = np.max(parameters["alpha_index_basic"]) + 1
    inv_dist_powers_ = np.zeros(max_alpha_index_basic)
    coords_powers_ = np.zeros( (max_alpha_index_basic, 3) )

    linear_mults = np.ones(parameters["alpha_scalar_moments"])
    max_linear = 1e-10 * np.ones(parameters["alpha_scalar_moments"])

    initializedVecs = {
                "max_dist":max_dist,
                "min_dist":mult, 
                "rb_vals":rb_vals,
                "rb_ders":rb_ders,
                "moment_vals":moment_vals,
                "basis_vals":basis_vals, 
                "site_energy_ders_wrt_moments_":site_energy_ders_wrt_moments_, 
                "max_alpha_index_basic":max_alpha_index_basic,
                "inv_dist_powers_":inv_dist_powers_, 
                "coords_powers_":coords_powers_,
                "linear_mults":linear_mults,
                "max_linear":max_linear,
                "mult":mult,
                }
    #
    return initializedVecs
#

def belongs(x, xmin, xmax):
    return (xmin <= x) and (x <= xmax)
#

# rb_vals, rb_ders = rb_Calc(r, min_dist, max_dist, mult, rb_vals, rb_ders, scaling, rb_size)
def rb_Calc(r, min_dist, max_dist, mult, rb_vals, rb_ders, scaling, rb_size):
    # from src/radial_basis.cpp: void RadialBasis_Chebyshev::RB_Calc(double r)
    
    # if parameters["rbasis_type"] == "RBChebyshev":
    assert belongs(r, min_dist, max_dist), "r does not belong to [min_dist, max_dist]: " + str(r) + " min_dist=" +str(min_dist) + " max_dist=" + str(max_dist)

    ksi = ( (2 * r) - (min_dist + max_dist)) / (max_dist - min_dist)
    R = r - max_dist
    R2 = R ** 2
	
    # scaling is not read in RadialBasis_Chebyshev::RB_Calc(), so it keeps its default value = 1
    rb_vals[0] = R2 # instead of `scaling * (1 * R2)``
	# rb_ders[0] = scaling * (0 * R2 + 2 * R) <<< why 0 * something??? I have to simplify below
    rb_ders[0] = 2.0 * R
    rb_vals[1] = ksi * R2
    rb_ders[1] = (mult * R2) + (2.0 * ksi * R)
    for i in range(2, rb_size):
        v1 = rb_vals[i - 1]
        v2 = rb_vals[i - 2]
        d1 = rb_ders[i - 1]
        d2 = rb_ders[i - 2]
        rb_vals[i] = (2.0 * ksi * v1) - v2
        rb_ders[i] = 2.0 * ( (mult * v1) + (ksi * d1) ) - d2
    #
    return rb_vals, rb_ders
#

# inv_dist_powers_, coords_powers_ = _pows(r, inv_dist_powers_, coords_powers_, max_alpha_index_basic, NeighbVect_j)
def _pows(r, inv_dist_powers_, coords_powers_, max_alpha_index_basic, NeighbVect_j):
    inv_dist_powers_[0] = 1.0
    coords_powers_[0]   = [1.0, 1.0, 1.0]
    for k in range(1, max_alpha_index_basic):
        inv_dist_powers_[k] = inv_dist_powers_[k - 1] / r
        for a in range(3):
            coords_powers_[k][a] = coords_powers_[k - 1][a] * NeighbVect_j[a]
    #    
    return inv_dist_powers_, coords_powers_
#

# buff_site_energy_ = calcSiteEnergyDers( nbh, type_central, # type of central atom at `i`
#                         linear_coeffs, lenNbh, 
#                         alpha_index_basic, alpha_index_basic_count, 
#                         max_alpha_index_basic,
#                         alpha_index_times, alpha_index_times_count,
#                         alpha_scalar_moments, alpha_moment_mapping,
#                         regression_coeffs, moment_vals,
#                         inv_dist_powers_, coords_powers_,
#                         species_count, radial_func_count, rb_size,
#                         min_dist, max_dist, mult, rb_vals, rb_ders, scaling,
#                         linear_mults, max_linear
#                         )
def calcSiteEnergyDers( nbh,
                        type_central, # type of central atom at `i`
                        linear_coeffs,
                        lenNbh, 
                        alpha_index_basic,
                        alpha_index_basic_count, 
                        max_alpha_index_basic,
                        alpha_index_times,
                        alpha_index_times_count,
                        alpha_scalar_moments,
                        alpha_moment_mapping,
                        regression_coeffs,
                        moment_vals,
                        inv_dist_powers_,
                        coords_powers_,
                        species_count,
                        radial_func_count,
                        rb_size,
                        min_dist,
                        max_dist,
                        mult,
                        rb_vals,
                        rb_ders,
                        scaling,
                        linear_mults,
                        max_linear
                        ):
    #
    # from dev_src/mtpr.cpp: void MLMTPR::CalcSiteEnergyDers(const Neighborhood& nbh)
    #

    buff_site_energy_ = 0.0

    # initialize
    moment_vals *= 0.0

    # lenNbh = len(nbh)

    # dicTypes = {"C":0, "O": 1}
    # types = [0,0,0,0, 1,1] #just an example

    # C = species_count   		#number of different species in current potential
    # K = radial_func_count		#number of radial functions in current potential
    # R = rb_size                 #number of Chebyshev polynomials constituting one radial function

    # moment_jacobian_ = np.zeros((alpha_index_basic_count, lenNbh, 3))

    assert type_central < species_count, "Too few species count in the MTP potential!"

    for j in range(lenNbh):
        # NeighbVect_j = nbh.vecs[j] #<<<<<<
        # r = nbh.dists[j]  #<<<<<<

        NeighbVect_j = nbh[j]["D"]
        r = nbh[j]["d"]

        rb_vals, rb_ders = rb_Calc(r, min_dist, max_dist, mult, rb_vals, rb_ders, scaling, rb_size)
        #

        rb_vals *= scaling # rb_vals is numpy array, so we can directly multyply by a float if array is float
        rb_ders *= scaling

        inv_dist_powers_, coords_powers_ = _pows(r, inv_dist_powers_, coords_powers_, max_alpha_index_basic, NeighbVect_j)

        # type_outer = nbh.types[j] #<<<<<<
        type_outer = nbh[j]["type"]

        for i in range(alpha_index_basic_count):
            val = 0.0
            # der = 0.0
            mu = alpha_index_basic[i, 0]

            # for l in range(rb_size):
            #     m = (type_central*C + type_outer)*K*R + mu * R + l
            #     val += regression_coeffs[m] * rb_vals[l]
            #     der += regression_coeffs[m] * rb_ders[l]

            for l in range(rb_size):
                coef = regression_coeffs[type_central, type_outer, mu, l]
                val += coef * rb_vals[l]
                # der += coef * rb_ders[l]
            #
            
            a1 = alpha_index_basic[i, 1]
            a2 = alpha_index_basic[i, 2]
            a3 = alpha_index_basic[i, 3]
            
            k = a1 + a2 + a3
            inv_powk = inv_dist_powers_[k]
            val *= inv_powk
            # der = (der * inv_powk) - (k * val / r)
            
            pow0 = coords_powers_[a1, 0]
            pow1 = coords_powers_[a2, 1]
            pow2 = coords_powers_[a3, 2]
            
            mult0 = pow0 * pow1 * pow2
            
            moment_vals[i] += val * mult0
            # mult0 *= der / r
            # moment_jacobian_[i, j, 0] += mult0 * NeighbVect_j[0]
            # moment_jacobian_[i, j, 1] += mult0 * NeighbVect_j[1]
            # moment_jacobian_[i, j, 2] += mult0 * NeighbVect_j[2]
            
            # if a1 != 0:
            #     prod = val * a1 * coords_powers_[a1 - 1, 0] * pow1 * pow2
            #     moment_jacobian_[i, j, 0] += prod
            # #
            # if a2 != 0:
            #     prod = val * a2 * pow0 * coords_powers_[a2 - 1, 1] * pow2
            #     moment_jacobian_[i, j, 1] += prod
			# #
            # if a3 != 0:
            #     prod = val * a3 * pow0 * pow1 * coords_powers_[a3 - 1, 2]
            #     moment_jacobian_[i, j, 2] += prod
			# #
        #

        ## Repulsive term
        ## I think it was not implemented in the C++ MTP code)
		## if (p_RadialBasis->GetRBTypeString() == "RBChebyshev_repuls")
        ## this seems arbitrary, I removed it:
        ## if r < min_dist:
		## 	multiplier = 10000;
		## 	buff_site_energy_ += multiplier*(exp(-10*(r-1)) - exp(-10*(min_dist-1)))
		## 	for (int a = 0; a < 3; a++)
		## 		buff_site_energy_ders_[j][a] += -10 * multiplier*(exp(-10 * (r - 1))/ nbh.dists[j])*nbh.vecs[j][a];
		# #
    #

    # Next: calculating non-elementary b_i
    for i in range(alpha_index_times_count):
        val0 = moment_vals[ alpha_index_times[i, 0] ] # float
        val1 = moment_vals[ alpha_index_times[i, 1] ] # float
        val2 = alpha_index_times[i, 2]                # integer
        moment_vals[alpha_index_times[i, 3]] += val2 * val0 * val1
	#

    # renewing maximum absolute values
    # for i in range(alpha_scalar_moments):
        # max_linear[i] = max(
        #     max_linear[i],
        #     abs(linear_coeffs[species_count + i] * moment_vals[alpha_moment_mapping[i]])
        #     )
        # #
    #

    # convolving with coefficients
    buff_site_energy_ += linear_coeffs[type_central]

    for i in range(alpha_scalar_moments):
        # buff_site_energy_ += linear_coeffs[species_count + i] * linear_mults[i] * moment_vals[alpha_moment_mapping[i]]
        # I simplified because linear_mults[i] = 1
        buff_site_energy_ += linear_coeffs[species_count + i] * moment_vals[alpha_moment_mapping[i]]
    #

    return buff_site_energy_
#

def CalcEFS(atoms, neighborhoods, type_centrals, params, vecs):
    # from src/basic_mlip.cpp:  void AnyLocalMLIP::CalcEFS(Configuration& cfg)
    energy = 0.0
    for i in range(len(atoms)):
        nbh = neighborhoods[i]
        lenNbh = len(nbh)
        type_central = type_centrals[i]

        # energy += calcSiteEnergyDers(nbh, type_central)
        # 
        energy += calcSiteEnergyDers(   nbh,
                                        type_central, # type of central atom at `i`
                                        params["linear_coeffs"],
                                        lenNbh,
                                        params["alpha_index_basic"],
                                        params["alpha_index_basic_count"],
                                        vecs["max_alpha_index_basic"],
                                        params["alpha_index_times"],
                                        params["alpha_index_times_count"],
                                        params["alpha_scalar_moments"],
                                        params["alpha_moment_mapping"],
                                        params["regression_coeffs"],
                                        vecs["moment_vals"],
                                        vecs["inv_dist_powers_"],
                                        vecs["coords_powers_"],
                                        params["species_count"],
                                        params["radial_func_count"],
                                        params["rb_size"],
                                        params["min_dist"],
                                        params["max_dist"],
                                        vecs["mult"],
                                        vecs["rb_vals"],
                                        vecs["rb_ders"],
                                        params["scaling"],
                                        vecs["linear_mults"],
                                        vecs["max_linear"]
                                    )        
        #
        # print(energy)
    #
    return energy
#


# neighborhoods = get_neighborhoods(atoms, cutOff)
def get_neighborhoods(atoms, listCutOffs, dictionaryTypes):
    neighborList = neighborlist.NeighborList(listCutOffs, self_interaction=False, bothways=True)
    neighborList.update(atoms)
    li, lj, ld, lD = neighborlist.neighbor_list('ijdD', atoms, listCutOffs)

    neighborhoods = {i:[] for i in range(len(atoms))}
    for k in range(len(li)):
        i = li[k]
        j = lj[k]
        d = ld[k]
        D = lD[k]
        t = dictionaryTypes[ atoms[j].symbol ]
        neighborhoods[i].append( { "j":j, "d":d, "D":D, "type":t } )
    #
    return neighborhoods
#

# type_centrals = get_type_centrals(atoms, dictionaryTypes)
def get_type_centrals(atoms, dictionaryTypes):
    # dictionaryTypes = { "C":0, "O":1 }
    type_centrals = [dictionaryTypes[i] for i in atoms.get_chemical_symbols()]
    return type_centrals
#