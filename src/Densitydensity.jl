ifexch(j::Int, k::Int, sys, range, exch::Float64) = ( 1 - (dis(j, k, sys; range=range) == 1) * exch )

""" A general principle is that we always use native indexing in the subsystem"""
DenDenNeighbor(sys::systems, j) = []

# function DenDenNeighbor(sys::QE_flat_SIAM, j)


#     _, chain_begin, chain_end = get_sys_loc(sys, j)
#     λ_ee, _, exch, _, range, _, ζ = CoulombParameters(sys)

#     ifexch(j, k, sys) = ( 1 - (dis(j, k, sys) == 1) * exch )

#     if chain_begin <= j <= chain_end

#         return [ [λ_ee * ifexch(j, k, sys) / ( dis(j, k, sys) + ζ), k] for k in max(chain_begin, j - range) : j - 1]


#     # center site interacting, since the sys is completely symmetric, the sites are the only thing thats being changed
#     elseif j == left(sys) + 1

#         denden = []

#         for k in 1:min(range, siteseach(sys))

#             # we calculate interact ONCE, as if we are at site 0
#             strength = center_ee(sys) * ifexch(0, k, sys) / (dis(0, k, sys) + ζ)
#             #left
#             for l in 1:legleft(sys)
#                 append!(denden, [[strength, j - k - (l -1) * (siteseach(sys) + QESITES)]])
#             end 

#             #right
#             for r in 1:legright(sys)
#                 append!(denden, [[strength, j + k + (r -1) * (siteseach(sys) + QESITES)]])
#             end 

#         end 

#         @show denden
#         return denden
#     else

#         return []
#     end 

# end 

# DenDenNeighbor(sys::QE_G_SIAM, j) = DenDenNeighbor(sys.system, j)

# function DenDenNeighbor(sys::QE_parallel, j)

#     uppertotal = get_uppertotal(sys)
#     systotal = get_systotal(sys)
#     center_lower = systotal - 1

#     if j <= uppertotal
#         den =  DenDenNeighbor(sys.upper, j)

#     elseif j < center_lower
#         den =  [ [U, k + uppertotal] for (U, k) in DenDenNeighbor(sys.lower, j - uppertotal)]

#     # the contact site, physically sits in the middle of two chains. We move to the native chain coordinates
#     elseif j == center_lower
#         upperchain = get_upperchain(sys)
#         lowerchain = get_lowerchain(sys)

#         uppercenter = div(upperchain, 2) + 0.5
#         lowercenter = div(lowerchain, 2) + 0.5
        
#         upperoffset = QESITES
#         loweroffset = uppertotal + QESITES

#         den = vcat(
#             [ [ center_ee(sys) / (dis(i, uppercenter, sys) + center_dis(sys)), i + upperoffset] for i in  max( 1, uppercenter - center_range(sys) ) : min( upperchain, uppercenter + center_range(sys))],
#             [ [ center_ee(sys) / (dis(i, lowercenter, sys) + center_dis(sys)), i + loweroffset] for i in  max( 1, lowercenter - center_range(sys) ) : min( lowerchain, lowercenter + center_range(sys))]
#         )

#     else

#         den = []
#     end 

#     return den
# end 

function DenDenNeighbor(sys::QE_two, j)

    systotal = get_systotal(sys)

    if j <3 || j > systotal - 2
        den = []

    else
        # shift the position of each j
        den = DenDenNeighbor(sys.chain_only, j - 2)

        # shift back position of each k
        den = [ [U, k + 2] for (U, k) in den]
    end

    return den
end 


function DenDenNeighbor(sys::QE_HOM, j)

    uppertotal = get_uppertotal(sys)

    minrange = max(QESITES + 1, ceil(true_center(sys) - center_range(sys)))
    maxrange = min(uppertotal - QESITES, floor(true_center(sys) + center_range(sys)))

    #@show minrange, maxrange

    if j <= uppertotal
        den =  DenDenNeighbor(sys.upper, j)

        # here, all the sites within the range in the upper chain are interacting with the lower chain wihtin that range window

        if minrange <= j <= maxrange
            coupling = [ [ center_ee(sys) / parallel_dis(i, j, sys) , i + uppertotal] for i in  minrange:maxrange]

            den = vcat(den, coupling)
        end 

    else
        den =  [ [U, k + uppertotal] for (U, k) in DenDenNeighbor(sys.lower, j - uppertotal)]
    end 

    return den
end 

function DenDenNeighbor(sys::Union{Chain_only, Rectangular}, j::Int; left_offset=0)

    λ_ee, _, exch, _, range, _, ζ = CoulombParameters(sys)

    adj_j = j - left_offset

    return [ [λ_ee * ifexch(adj_j, k, sys, range, exch) / ( dis(adj_j, k, sys; range=range) + ζ), k + left_offset] for k in 1:adj_j- 1]
    
end 

DenDenNeighbor(sys::biased_chain, j::Int, left_offset=0) = sys.chain_start <= j - left_offset < sys.chain_start + L(sys.chain) ? DenDenNeighbor(sys.chain, j, left_offset=left_offset + (sys.chain_start - 1)) : []

DenDenNeighbor(sys::GQS, j::Int) = DenDenNeighbor(sys.chain_only, j)

function DenDenNeighbor(sys::DPT, j::Int)

    # if j == L(sys) + 1
    if j == dd_lower(sys)

        #den =  vcat( [ [U(sys), j - k] for k in 1:couple_range(sys)], [[U(sys), j + 1 + k] for k in 1:couple_range(sys)])

        den =  vcat( [ [U(sys), k] for k in L_contact(sys):L_end(sys)], [[U(sys), k] for k in R_begin(sys):R_contact(sys)])
    else
        den = []
    end 

    return den

end 


function DenDenNeighbor(sys::DPT_mixed, j::Int)

    # if center region is not under mixed basis, they are still interacting spatially 
    if !includeU(sys)
        den = DenDenNeighbor(sys.dpt, j)
    else
        den = []
    end 

    return den
end 

DenDenNeighbor(res::reservoir, j::Int; left_offset=0) = []

function DenDenNeighbor(sys::SD_array, j::Int)

    source = get_systotal(sys.source)
    array = get_systotal(sys.array)

    if j <= source
        return DenDenNeighbor(sys.source, j)

    elseif j <= source + array
        return DenDenNeighbor(sys.array, j; left_offset=source)

    else
        return DenDenNeighbor(sys.drain, j; left_offset= source + array)

    end 

end 


function DenDenNeighbor(sys::DPT_avg, j::Int)

    DenDenNeighbor(sys.dpt, j)

end 


DenDenNeighbor(sys::DPT_graph, j::Int) = DenDenNeighbor(sys.dpt, j)

ddoperators(sys::systems) = systype(sys) == "Fermion" ? [["N", "N"]] : [["Ntot", "Ntot"]]
ddoperators(sys::DPT_avg) =  [["Ndn", "Nup"]]

function add_DensityDensity!(sys::systems, res::OpSum)
    

    @info "Adding DemDen"

    for j=1 :get_systotal(sys)
        # E-E and N-E

        for (U..., k) in DenDenNeighbor(sys, j)
            
            k = trunc(Int, k)
            for (i, operators) in enumerate( ddoperators(sys))
                
                op1, op2 = operators
                res += U[i], op1, sitemap(sys, j), op2, sitemap(sys, k)

            end 

        end

    end 

    return res

end 

