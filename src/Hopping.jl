HoppingNeighbor(sys::systems, j::Int; left_offset=0) = []


# for a flat chain, return next NN until end
HoppingNeighbor(sys::Chain_only, j::Int; left_offset=0) = 0 < (j - left_offset) < L(sys) ? [[t(sys)..., j + 1] ] : []

HoppingNeighbor(sys::biased_chain, j::Int; left_offset=0) = sys.chain_start <= j - left_offset < sys.chain_start + L(sys.chain) ? HoppingNeighbor(sys.chain, j, left_offset=left_offset + (sys.chain_start - 1)) : []

HoppingNeighbor(sys::GQS, j ::Int) = HoppingNeighbor(sys.chain_only, j)
"""
Here, we have either the left_offset hopping in the QE, or the chain hopping
"""
function HoppingNeighbor(sys::QE_two, j::Int; left_offset=0)

    adj_j = j - left_offset
    # we only hop onwards, no hop in QE and last site
    if adj_j  <3 || adj_j > get_systotal(sys) - 3
        hop = []
    else
        hop = [[t(sys)..., j + 1]]
    end 

    return hop

end 

# function HoppingNeighbor(sys::QE_parallel, j::Int)

#     uppertotal = get_uppertotal(sys)
#     systotal = get_systotal(sys)

#     center_lower = systotal - 1

#     if j <= uppertotal
#         hop = HoppingNeighbor(sys.upper, j)

#     elseif j < center_lower
#         hop = [ [t, site + uppertotal] for (t, site) in HoppingNeighbor(sys.lower, j - uppertotal)]

#     elseif j == center_lower
#         hop = [ [center_internal_t(sys), systotal]]

#     else
#         hop = []
#     end 

#     return hop
# end 

# """Flattened X QE, each 'arm' is staggered, as QE + chain, ... , center, chain, QE, ...."""
# function HoppingNeighbor(sys::QE_flat_SIAM, j::Int)

#     hop = []
#     # determine the total number of sites in left
#     @show qe_loc, chain_begin, chain_end = get_sys_loc(sys, j)

#     # hop to next
#     if chain_begin <= j < chain_end
#         append!(hop, [[t(sys), j + 1]])
#     end 

#     # left end, to center
#     # right begin, to center
#     if (j == chain_end && j > qe_loc)  ||  (j == chain_begin && j < qe_loc)
#         append!(hop, [[center_t(sys), left(sys) + 1]])
#     end 

#     return hop


# end 

# HoppingNeighbor(sys::QE_G_SIAM, j::Int) = HoppingNeighbor(sys.system, j::Int)

function HoppingNeighbor(sys::QE_HOM, j::Int)

    uppertotal = get_uppertotal(sys)
    systotal = get_systotal(sys)


    if j <= uppertotal
        hop = HoppingNeighbor(sys.upper, j)

    elseif j < systotal
        hop = HoppingNeighbor(sys.lower, j; left_offset=uppertotal)

    else
        hop = []
    end 

    return hop
end 



function HoppingNeighbor(sys::NF_square, j::Int)

    hop = []

    # not at end of col
    if j % L(sys) != 0
        append!(hop, [[t(sys)..., j + 1]])
    end 

    # not at end of row
    if div(j - 1, L(sys)) + 1 < L(sys)
        append!(hop, [[t(sys)..., j + L(sys)]])
    end 


    return hop

end 

function HoppingNeighbor(sys::Rectangular, j::Int; left_offset=0)

    hop = []
    adj_j = j - left_offset

    # not at end of col
    if adj_j % Lx(sys) != 0
        append!(hop, [[t(sys)..., j + 1]])
    end 

    # not at end of row
    if adj_j <= get_systotal(sys) - Lx(sys)
        append!(hop, [[t(sys)..., j + Lx(sys) ]])
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

    if L_begin(sys) <= j < L_end(sys) || R_begin(sys) <= j < R_end(sys)
        hop = [[t_reservoir(sys), j + 1]]

    elseif j == L_end(sys)
        hop =  [[t_reservoir(sys), R_begin(sys)]]

    elseif j == dd_lower(sys)
        hop = [[t_doubledot(sys), j + 1]]

    else
        hop = []
    end 

    return hop

end 

# no hopping except DD and center region (if applicable)
function HoppingNeighbor(sys::DPT_mixed, j::Int)

    if j == dd_lower(sys)
        hop = [[t_doubledot(sys), j + 1]]

    else
        # we connect within the spatial region
        if !includeU(sys) 
            
            # we connect within L
            if L_end(sys) - couple_range(sys) < j < L_end(sys)
                hop = [[t_reservoir(sys), j + 1]]
            
            # we connect within R
            elseif R_begin(sys) <= j < R_begin(sys) + couple_range(sys) - 1
                hop  = [[t_reservoir(sys), j + 1]]

            # we connect LR
            elseif j == L_end(sys)
                hop =  [[t_reservoir(sys), R_begin(sys)]]

            else
                hop = []
            end 

        else
            hop = []
        end 
    end 

    

    return hop

end 

HoppingNeighbor(sys::DPT_graph, j::Int) = HoppingNeighbor(sys.dpt, j)

function HoppingNeighbor(sys::DPT_avg, j::Int)

    if j == dd_lower(sys)
        return [[t_reservoir(sys), t_doubledot(sys), j + 1]]

    else
        hop =  HoppingNeighbor(sys.dpt, j)

        if !isempty(hop)
            return [[h[1], 0, h[2]] for h in hop]
        else
            return hop
        end 


    end 
    

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



function HoppingNeighbor(res::reservoir_spatial, j::Int, contacts::Array; left_offset=0)

    hop = []
    adj_j = j - left_offset

    # check if hopping to sys
    if adj_j == res.contact
        append!(hop, contacts)
    end 

    # check if can hop to NN
    if adj_j < get_systotal(res) 
        append!(hop, [[res.t..., j + 1 ]] )
    end 

    return hop

end 


function HoppingNeighbor(sys::SD_array, j::Int)

    source = get_systotal(sys.source)
    array = get_systotal(sys.array)

    if j <= source
        return HoppingNeighbor(sys.source, j, sys.s_contacts)

    elseif j <= source + array
        return HoppingNeighbor(sys.array, j; left_offset=source)

    else
        return HoppingNeighbor(sys.drain, j, sys.d_contacts; left_offset= source + array)

    end 
    
end 



"""This function only concerns with the 'simple' hopping, complex terms such as QE dipole offsets are calculated elsewhere"""
function add_hop!(sys::systems, res::OpSum)
    

    @info "Adding all hopping"
    sys_type = systype(sys)
    systotal = get_systotal(sys)

    if sys_type == "Fermion"
        operators = [ ["C", "Cdag"]]

    elseif sys_type == "Electron"
        operators = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]

    end 

    for j in 1:systotal

        for (v..., k) in HoppingNeighbor(sys, j)

            k = trunc(Int, k)
            for (i, operator) in enumerate(operators)

                if v[i] != 0
                    op1, op2 = operator
                    

                    res += v[i], op1, sitemap(sys, j), op2, sitemap(sys, k)
                    res += v[i], op1, sitemap(sys, k), op2, sitemap(sys, j)
                end 
            end 

        end 
    end 
    
    return res

end 
