# for a flat chain, return next NN until end
HoppingNeighbor(sys::Chain_only, j::Int) = j < L(sys) ? [[t(sys), j + 1] ] : []


"""
Here, we have either the offset hopping in the QE, or the chain hopping
"""
function HoppingNeighbor(sys::QE_two, j::Int)


    # we only hop onwards, no hop in QE and last site
    if j <3 || j > get_systotal(sys) - 3
        hop = []
    else
        hop = [[t(sys), j + 1]]
    end 

    return hop

end 

function HoppingNeighbor(sys::QE_parallel, j::Int)

    uppertotal = get_uppertotal(sys)

    if j <= uppertotal
        hop = HoppingNeighbor(sys.upper, j)

    elseif j < get_systotal(sys)
        hop = [ [t, site + uppertotal] for (t, site) in HoppingNeighbor(sys.lower, j - uppertotal)]

    else
        hop = []
    end 

    return hop
end 

"""Flattened X QE, each 'arm' is staggered, as QE + chain, ... , center, chain, QE, ...."""
function HoppingNeighbor(sys::QE_flat_SIAM, j::Int)

    hop = []
    # determine the total number of sites in left
    @show qe_loc, chain_begin, chain_end = get_sys_loc(sys, j)

    # hop to next
    if chain_begin <= j < chain_end
        append!(hop, [[t(sys), j + 1]])
    end 

    # left end, to center
    # right begin, to center
    if (j == chain_end && j > qe_loc)  ||  (j == chain_begin && j < qe_loc)
        append!(hop, [[center_t(sys), left(sys) + 1]])
    end 

    return hop


end 

HoppingNeighbor(sys::QE_G_SIAM, j::Int) = HoppingNeighbor(sys.system, j::Int)

function HoppingNeighbor(sys::NF_square, j::Int)

    hop = []

    # not at end of col
    if j % L(sys) != 0
        append!(hop, [[t(sys), j + 1]])
    end 

    # not at end of row
    if div(j - 1, L(sys)) + 1 < L(sys)
        append!(hop, [[t(sys), j + L(sys)]])
    end 


    return hop

end 

function HoppingNeighbor(sys::DPT, j::Int)

    # # if L or R, no contact
    # if j < L(sys)  || (j > L(sys) + 2 && j < get_systotal(sys))
    #     hop = [[t_reservoir(sys), j + 1]]

    # # contact
    # elseif j == L(sys) && contact(sys)
    #     hop = [[contact_t(sys), j + 3]]

    # # lower dot
    # elseif j == L(sys) + 1
    #     hop = [[t_doubledot(sys), j + 1]]

    if j < L(sys) + R(sys)
        hop = [[t_reservoir(sys), j + 1]]

    elseif j == L(sys) + R(sys) + 1
        hop = [[t_doubledot(sys), j + 1]]

    else
        hop = []
    end 

    return hop

end 

function HoppingNeighbor(sys::LSR_SIAM, j::Int)

    # if L or R
    if j < L(sys)  || (j > L(sys) + 1 && j < get_systotal(sys))
        hop = [[t_reservoir(sys), j + 1]]

    elseif j == L(sys) || j == L(sys) + 1
        hop = [[t_couple(sys), j + 1]]

    else
        hop = []
    end 


    return hop
end 

"""This function only concerns with the 'simple' hopping, complex terms such as QE dipole offsets are calculated elsewhere"""
function add_hop!(sys::systems, res::OpSum)
    

    println("Adding all hopping")
    sys_type = type(sys)
    systotal = get_systotal(sys)

    if sys_type == "Fermion"
        operators = [ ["C", "Cdag"]]

    elseif sys_type == "Electron"
        operators = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]

    end 

    for j in 1:systotal

        for (v, k) in HoppingNeighbor(sys, j)

            k = trunc(Int, k)
            for operator in operators

                op1, op2 = operator

                res += v, op1, sitemap(sys, j), op2, sitemap(sys, k)
                res += v, op1, sitemap(sys, k), op2, sitemap(sys, j)
            end 

        end 
    end 
    
    return res

end 
