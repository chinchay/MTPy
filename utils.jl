#%%
function getParam(iLine)
    return strip( split(lines[iLine], "=")[2] , ' ')
end

function getParam2(n, iLine, typeNumber)
    v = zeros(n)
    l = split( strip( getParam(iLine), ['{', '}'] ), ',' )
    for i in 1:n
        if typeNumber == "Int"
            v[i] = parse(Int, l[i])
        elseif typeNumber == "Float64"
            v[i] = parse(Float64, l[i])
        end
    end
    return v
end    

function getParam3(n, iLine)
    v = zeros(Int, (n, 4))
    l = split( strip( getParam(iLine), ['{', '}'] ), "}, {"  )
    for i in 1:n
            temp = split( l[i], ',' )
            for j in 1:4
                v[i, j] = parse(Int, temp[j])
            end
    end
    return v
end    

function load(filename="workdir/pot.mtp")
    lines = readlines(filename)

    @assert lines[1] == "MTP" "Can read only MTP format potentials"

    version = getParam(2)
    @assert version == "1.1.0" "MTP file must have version \"1.1.0\""

    potential_name     = getParam(3)
    scaling            = parse(Float64, getParam(4))
    species_count      = parse(Int, getParam(5))
    potential_tag      = getParam(6)
    radial_basis_type  = getParam(7)
    min_dist           = parse(Float64, getParam(8))
    max_dist           = parse(Float64, getParam(9))
    radial_basis_size  = parse(Int, getParam(10))
    radial_funcs_count = parse(Int, getParam(11))
    
    iline = 13
    regression_coeffs = zeros( (species_count, species_count, radial_funcs_count, radial_basis_size) )
    for s1 in 1:species_count
        for s2 in 1:species_count
            iline += 1
            for i in 1:radial_funcs_count
                list_t = split( strip( lines[15], ['\t', '{', '}'] ), ',' )
                iline += 1
                for j in 1:radial_basis_size
                    regression_coeffs[s1, s2, i, j] = parse(Float64, list_t[j])
                end
            end
        end
    end
    #
    
    alpha_moments_count = parse(Int, getParam(iline))
    
    alpha_index_basic_count = parse(Int, getParam(iline + 1))
    alpha_index_basic       = getParam3(alpha_index_basic_count, iline + 2)
    
    radial_func_max = maximum([ alpha_index_basic[i, 1] for i in 1:alpha_index_basic_count ])
    @assert radial_func_max == (radial_funcs_count - 1) "Wrong number of radial functions specified"

    alpha_index_times_count = parse(Int, getParam(iline + 3))
    alpha_index_times       = getParam3(alpha_index_times_count, iline + 4)
    
    alpha_scalar_moments = parse(Int, getParam(iline + 5))
    alpha_moment_mapping = getParam2(alpha_scalar_moments, iline + 6, "Int")

    alpha_count = alpha_scalar_moments + 1
    
    species_coeffs       = getParam2(species_count, iline + 7, "Float64")
    moment_coeffs        = getParam2(alpha_scalar_moments, iline + 8, "Float64")
    #
    params = Dict(
            "scaling" => scaling,
            "species_count" => species_count,
            "rbasis_type" => radial_basis_type,
            "min_dist" => min_dist,
            "max_dist" => max_dist,
            "rb_size" => radial_basis_size,
            "radial_func_count" => radial_funcs_count,
            "regression_coeffs" => regression_coeffs,
            "alpha_moments_count" => alpha_moments_count,
            "alpha_index_basic_count" => alpha_index_basic_count,
            "alpha_index_basic" => alpha_index_basic,
            "alpha_index_times" => alpha_index_times,
            "alpha_scalar_moments" => alpha_scalar_moments,
            "alpha_count" => alpha_count,
            # "linear_coeffs":linear_coeffs,
            "species_coeffs" => species_coeffs,
            "moment_coeffs" => moment_coeffs,
            "alpha_index_times_count" => alpha_index_times_count,
            "alpha_moment_mapping" => alpha_moment_mapping,
            )
    #
    return params
end

function init_vecs(parameters)
    rb_vals     = zeros(parameters["rb_size"])
    rb_ders     = zeros(parameters["rb_size"])
    moment_vals = zeros(parameters["alpha_moments_count"])
    basis_vals  = zeros(parameters["alpha_count"])
    site_energy_ders_wrt_moments_ = zeros(parameters["alpha_moments_count"])
    
    mult        = 2.0 / (parameters["max_dist"] - parameters["min_dist"])
    
    max_alpha_index_basic = maximum(parameters["alpha_index_basic"]) + 1
    inv_dist_powers_ = zeros(max_alpha_index_basic)
    coords_powers_   = zeros( (max_alpha_index_basic, 3) )

    linear_mults = ones(parameters["alpha_scalar_moments"])
    max_linear = 1e-10 * ones(parameters["alpha_scalar_moments"])

    # # calculate moment_vals[i] for i only in the following list:
    # alpha_index_basic_count = parameters["alpha_index_basic_count"]
    # alpha_index_times_count = parameters["alpha_index_times_count"]
    # alpha_index_times = parameters["alpha_index_times"]
    # alpha_moment_mapping = parameters["alpha_moment_mapping"]
    # l1 = [ i for i in range(alpha_index_basic_count) if i in alpha_moment_mapping ]
    # l2 = [ alpha_index_times[i, 0] for i in range(alpha_index_times_count) if alpha_index_times[i, 0] < alpha_index_basic_count ]
    # l3 = [ alpha_index_times[i, 1] for i in range(alpha_index_times_count) if alpha_index_times[i, 1] < alpha_index_basic_count ]
    # l4 = [ alpha_index_times[i, 3] for i in range(alpha_index_times_count) if alpha_index_times[i, 3] < alpha_index_basic_count ]    
    # l1 = l1 + list( set(l2) - set(l1) )
    # l1 = l1 + list( set(l3) - set(l1) )
    # i_moment_vals = l1 + list( set(l4) - set(l1) )

    initializedVecs = Dict(
                "rb_vals" => rb_vals,
                "rb_ders" => rb_ders,
                "moment_vals" => moment_vals,
                "basis_vals" => basis_vals, 
                "site_energy_ders_wrt_moments_" => site_energy_ders_wrt_moments_, 
                "max_alpha_index_basic" => max_alpha_index_basic,
                "inv_dist_powers_" => inv_dist_powers_, 
                "coords_powers_" => coords_powers_,
                "linear_mults" => linear_mults,
                "max_linear" => max_linear,
                "mult" => mult
    )
    #
    return initializedVecs
end

function correctIndexesInParameters!(parameters)
    # mutates `parameters`
    # correcting to index starting in 1 in Julia, instead of 0 as in C++ (MTP code)
    for i in 1:parameters["alpha_index_basic_count"]
        parameters["alpha_index_basic"][i, 1] += 1 # = mu
    end
    for i in 1:parameters["alpha_index_times_count"]
        parameters["alpha_index_times"][i, 1] += 1 # val0 = moment_vals[ alpha_index_times[i, 0] ] in finish_moment_vals()
        parameters["alpha_index_times"][i, 2] += 1 # val1 = moment_vals[ alpha_index_times[i, 1] ]
        parameters["alpha_index_times"][i, 4] += 1 # moment_vals[alpha_index_times[i, 3]] += val2 * val0 * val1
    end
    for i in 1:parameters["alpha_scalar_moments"]
        parameters["alpha_moment_mapping"][i] += 1
    end
    parameters["alpha_index_basic"] = alpha_index_basic
    parameters["alpha_index_times"] = alpha_index_times
    parameters["alpha_moment_mapping"] = alpha_moment_mapping
end

function belongs(x, xmin, xmax)
    return (xmin <= x) && (x <= xmax)
end

# rb_Calc!(r, min_dist, max_dist, mult, rb_vals, rb_ders, scaling, rb_size)
function rb_Calc!(r, min_dist, max_dist, mult, rb_vals, rb_ders, scaling, rb_size)
    # mutates rb_vals, rb_ders
    # from src/radial_basis.cpp: void RadialBasis_Chebyshev::RB_Calc(double r)
    
    # if parameters["rbasis_type"] == "RBChebyshev":
    @assert belongs(r, min_dist, max_dist) "r does not belong to 
                [min_dist, max_dist]: " + string(r) + " min_dist=" + 
                string(min_dist) + " max_dist=" + string(max_dist)
    #

    ksi = ( (2 * r) - (min_dist + max_dist)) / (max_dist - min_dist)
    R   = r - max_dist
    R2  = R ^ 2
	
    # scaling is not read in RadialBasis_Chebyshev::RB_Calc(), so it keeps its default value = 1
    rb_vals[1] = R2 # instead of `scaling * (1 * R2)``
	# rb_ders[0] = scaling * (0 * R2 + 2 * R) <<< why 0 * something??? I have to simplify below
    rb_ders[1] = 2.0 * R

    rb_vals[2] = ksi * R2
    rb_ders[2] = (mult * R2) + (2.0 * ksi * R)

    for i in 3:rb_size
        v1 = rb_vals[i - 1]
        v2 = rb_vals[i - 2]
        d1 = rb_ders[i - 1]
        d2 = rb_ders[i - 2]
        rb_vals[i] = (2.0 * ksi * v1) - v2
        rb_ders[i] = 2.0 * ( (mult * v1) + (ksi * d1) ) - d2
    end
    # return rb_vals, rb_ders
end

# _pows(r, inv_dist_powers_, coords_powers_, max_alpha_index_basic, NeighbVect_j)
function _pows!(r, inv_dist_powers_, coords_powers_, max_alpha_index_basic, NeighbVect_j)
    # mutates inv_dist_powers_, coords_powers_
    inv_dist_powers_[1] = 1.0
    coords_powers_[1]   = [1.0, 1.0, 1.0]
    for k in 2:max_alpha_index_basic
        inv_dist_powers_[k] = inv_dist_powers_[k - 1] / r
        for i in 1:3
            coords_powers_[k][i] = coords_powers_[k - 1][i] * NeighbVect_j[i]
        end
    end
    return inv_dist_powers_, coords_powers_
end

# Next: calculating non-elementary b_i
function finish_moment_vals!(
                moment_vals,
                alpha_index_times_count,
                alpha_index_times,
                )
    # mutates moment_vals
    for i in 1:alpha_index_times_count
        val0 = moment_vals[ alpha_index_times[i, 1] ] # float
        val1 = moment_vals[ alpha_index_times[i, 2] ] # float
        val2 = alpha_index_times[i, 3]                # integer
        moment_vals[alpha_index_times[i, 4]] += val2 * val0 * val1
	end
end

# val, der = get_val_der( mu, regression_coeffs, type_central, type_outer, rb_vals, rb_ders )
function get_val_der(    mu,
                    regression_coeffs,
                    type_central,
                    type_outer,
                    rb_vals,
                    rb_ders
                )
    #
    val = regression_coeffs[type_central, type_outer, mu, :] .* rb_vals
    der = regression_coeffs[type_central, type_outer, mu, :] .* rb_ders
    return val, der
end

function update_moment_vals!(
                    i_moment_vals,
                    moment_vals,
                    alpha_index_basic_count,
                    alpha_index_basic,
                    inv_dist_powers_,
                    regression_coeffs,
                    type_central,
                    type_outer,
                    rb_vals,
                    rb_ders,
                    coords_powers_
                    )
    # update_moment_vals
    #
    for i in 1:alpha_index_basic_count
    # for i in i_moment_vals:
        # print(i)
        val = 0.0
        # der = 0.0
        mu = alpha_index_basic[i, 1]

        # for l in range(rb_size):
        #     m = (type_central*C + type_outer)*K*R + mu * R + l
        #     val += regression_coeffs[m] * rb_vals[l]
        #     der += regression_coeffs[m] * rb_ders[l]

        # for l in range(rb_size):
        #     coef = regression_coeffs[type_central, type_outer, mu, l]
        #     val += coef * rb_vals[l]
        #     # der += coef * rb_ders[l]
        # #
        # depend on `r`:
        val, der = get_val_der( mu, regression_coeffs, 
                                type_central, type_outer, 
                                rb_vals, rb_ders )
        #
        
        a1 = alpha_index_basic[i, 2]
        a2 = alpha_index_basic[i, 3]
        a3 = alpha_index_basic[i, 4]
        
        k = a1 + a2 + a3
        inv_powk = inv_dist_powers_[k]
        val *= inv_powk
        # der = (der * inv_powk) - (k * val / r)
        
        pow0 = coords_powers_[a1, 1]
        pow1 = coords_powers_[a2, 2]
        pow2 = coords_powers_[a3, 3]
        
        mult0 = pow0 * pow1 * pow2
        
        
        
        #
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
    end
end

function calcSiteEnergyDers(
                        nbh,
                        type_central, # type of central atom at `i`
                        species_coeffs,
                        moment_coeffs,
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
                        max_linear,
                        i_moment_vals
                        )
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

    @assert type_central < species_count "Too few species count in the MTP potential!"

    for j in 1:lenNbh
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

        moment_vals = update_moment_vals(
                            i_moment_vals,
                            moment_vals,
                            alpha_index_basic_count,
                            alpha_index_basic,
                            inv_dist_powers_,
                            regression_coeffs,
                            type_central,
                            type_outer,
                            rb_vals,
                            rb_ders,
                            coords_powers_
                            )
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
    end

    # # Next: calculating non-elementary b_i
    # for i in range(alpha_index_times_count):
    #     val0 = moment_vals[ alpha_index_times[i, 0] ] # float
    #     val1 = moment_vals[ alpha_index_times[i, 1] ] # float
    #     val2 = alpha_index_times[i, 2]                # integer
    #     moment_vals[alpha_index_times[i, 3]] += val2 * val0 * val1
	# #
    moment_vals = finish_moment_vals(moment_vals, alpha_index_times_count, alpha_index_times)

    # renewing maximum absolute values
    # for i in range(alpha_scalar_moments):
        # max_linear[i] = max(
        #     max_linear[i],
        #     abs(linear_coeffs[species_count + i] * moment_vals[alpha_moment_mapping[i]])
        #     )
        # #
    #

    # convolving with coefficients
    buff_site_energy_ += species_coeffs[type_central]

    for i in 1:alpha_scalar_moments
        # buff_site_energy_ += linear_coeffs[species_count + i] * linear_mults[i] * moment_vals[alpha_moment_mapping[i]]
        # I simplified because linear_mults[i] = 1
        buff_site_energy_ += moment_coeffs[i] * moment_vals[alpha_moment_mapping[i]]
    end

    return buff_site_energy_
end

function CalcEFS(atoms, neighborhoods, type_centrals, params, vecs)
    # from src/basic_mlip.cpp:  void AnyLocalMLIP::CalcEFS(Configuration& cfg)
    energy = 0.0
    for i in 1:atoms.get_global_number_of_atoms()
        nbh = neighborhoods[i]
        lenNbh = len(nbh)
        type_central = type_centrals[i]

        # energy += calcSiteEnergyDers(nbh, type_central)
        # 
        energy += calcSiteEnergyDers(   nbh,
                                        type_central, # type of central atom at `i`
                                        params["species_coeffs"],
                                        params["moment_coeffs"],
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
                                        vecs["max_linear"],
                                        vecs["i_moment_vals"],
                                    )        
        #
        # print(energy)
    end
    return energy
end

function get_neighborhoods(atoms, listCutOffs, dictionaryTypes)
    neighborList = neighborlist.NeighborList(listCutOffs, self_interaction=False, bothways=True)
    neighborList.update(atoms)
    li, lj, ld, lD = neighborlist.neighbor_list("ijdD", atoms, listCutOffs)

    neighborhoods = Dict( i => [] for i in range(len(atoms)) )
    # neighborhoods = [ [] for i in range(len(atoms)) ]
    for k in 1:len(li)
        i = li[k]
        j = lj[k]
        d = ld[k]
        D = lD[k]
        t = dictionaryTypes[ atoms[j].symbol ]
        neighborhoods[i].append( Dict( "j"=>j, "d"=>d, "D"=>D, "type"=>t ) )
        # neighborhoods[i].append( [ j, d, D, t ] )
    end
    return neighborhoods
end

# type_centrals = get_type_centrals(atoms, dictionaryTypes)
function get_type_centrals(atoms, dictionaryTypes)
    # dictionaryTypes = { "C":0, "O":1 }
    type_centrals = [dictionaryTypes[i] for i in atoms.get_chemical_symbols()]
    return type_centrals
end

#%%

# x = 3

#%%
using ASE
atoms = bulk("Cu", cubic=true) * 2        # generate periodic Cu supercell
# deleteat!(at, 1)                       # vacancy defect
# emt = pyimport("ase.calculators.emt")  # import the EMT model
# calc = ASECalculator(emt.EMT())        # wrap it into a Julia Object
# @show energy(calc, at)

# atoms.get_global_number_of_atoms()
# len(atoms)

#%%
using JuLIP
at = bulk(:Si, cubic=true)



