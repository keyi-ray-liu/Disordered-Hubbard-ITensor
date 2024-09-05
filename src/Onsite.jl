Onsite(sys::systems, j; left_offset=0) = 0.0


function Onsite(sys::QE_two, j; left_offset=0) :: Float64

    adj_j = j - left_offset
    systotal = get_systotal(sys)

    if adj_j == 2 || adj_j == systotal
        onsite = 0.0

    elseif adj_j == 1 || adj_j == systotal - 1
        onsite = QEen(sys)

    else
        # adjust position of j
        onsite = Onsite(sys.chain, j; left_offset=2)
        
        # see if we have confining potential
        start = confine_start(sys)
        range = confine_range(sys)
        potential = confine_potential(sys)

        if 2 + start <= adj_j < 2 + start + range
            onsite += potential
        end 
    end

    return onsite
end 

# function Onsite(sys::QE_parallel, j) 

#     uppertotal = get_uppertotal(sys)
#     systotal = get_systotal(sys)
#     center_lower = systotal - 1

#     if j <= uppertotal
#         onsite = Onsite(sys.upper, j)

#     elseif j < center_lower
#         onsite = Onsite(sys.lower, j - uppertotal)

#     elseif j == center_lower
#         upperchain = get_upperchain(sys)
#         lowerchain = get_lowerchain(sys)

#         uppercenter = div(upperchain, 2) + 0.5
#         lowercenter = div(lowerchain, 2) + 0.5

#         onsite =  - ( sum(
#             [ center_ne(sys) / (dis(i, uppercenter, sys) + center_dis(sys)) for i in  max( 1, uppercenter - center_range(sys) ) : min( upperchain, uppercenter + center_range(sys))]) + 
#             sum( [ center_ne(sys) / (dis(i, lowercenter, sys) + center_dis(sys)) for i in  max( 1, lowercenter - center_range(sys) ) : min( lowerchain, lowercenter + center_range(sys))]
#         ) ) 

#     else
#         onsite = 0.0
#     end 

#     return onsite
    
# end 

# function Onsite(sys::QE_flat_SIAM, chain_begin, chain_end, j)

#     _, λ_ne, _, _, range, CN, ζ = CoulombParameters(sys)

#     λ_ne *= CN

#     if j == left(sys) + 1

#         onsite =  sum([ - center_ne(sys) * CN / (dis(0, k, sys) + ζ) for k in 1: min(range, siteseach(sys))]) * (legleft(sys) + legright(sys))


#     else
        
#         onsite = - λ_ne * CN *  ( sum([ 1/(dis(j, k, sys) + ζ) for k in max(chain_begin, j - range) : min(chain_end, j + range)])   - 1/ζ)
#     end 

#     return onsite

# end 


# function Onsite(sys::QE_flat_SIAM, j) :: Float64

#     qe_loc, chain_begin, chain_end = get_sys_loc(sys, j)

#     if j == qe_loc + 1 
#         onsite = 0.0

#     elseif j == qe_loc
#         onsite = QEen(sys)
#     else
#         # adjust position of j
#         onsite = Onsite(sys, chain_begin, chain_end, j)
#     end

#     return onsite
# end 

# Onsite(sys::QE_G_SIAM, j) = Onsite(sys.system, j)


function Onsite(sys::QE_HOM, j) 

    uppertotal = get_uppertotal(sys)
    minrange = max(QESITES + 1, ceil(true_center(sys) - center_range(sys)))
    maxrange = min(uppertotal - QESITES, floor(true_center(sys) + center_range(sys)))

    #@show minrange, maxrange

    if j <= uppertotal
        onsite = Onsite(sys.upper, j)

        if minrange <= j <= maxrange
            onsite -= sum([  center_ne(sys) / parallel_dis(i, j, sys) for i in  minrange:maxrange])
        end 

    else
        onsite = Onsite(sys.lower,j; left_offset=uppertotal)

        if minrange <= j - uppertotal <= maxrange
            onsite -= sum([  center_ne(sys) / parallel_dis(i, j - uppertotal, sys) for i in  minrange:maxrange])
        end 
    end 

    return onsite
    
end 



function Onsite(sys::Union{Rectangular, Chain}, j::Int; left_offset=0) :: Float64

    _, λ_ne, _, _, range, CN, ζ = CoulombParameters(sys)

    λ_ne *= CN
    systotal = get_systotal(sys)
    
    adj_j = j - left_offset
    onsite = - λ_ne * ( sum([ 1/(dis(adj_j, k, sys; range=range) + ζ) for k in 1:systotal])   - 1/ζ)

    return onsite

end 

# we biase everywhere else on the chain with a large positive potential
Onsite(sys::biased_chain, j::Int, left_offset=0) = sys.chain_start <= j - left_offset < sys.chain_start + L(sys.chain) ? Onsite(sys.chain, j, left_offset=left_offset + (sys.chain_start - 1)) : 500.0


Onsite(sys::SSH_chain, j::Int; left_offset =0) = Onsite(sys.chain, j; left_offset=left_offset)

Onsite(sys::GQS, j::Int) = Onsite(sys.chain, j)

"""For NF NxN, we add to the inner square, that is from row 2 to row N - 1, from col 2 to col N - 1"""
function Onsite(sys::NF_square, j::Int) :: Float64
    
    row = div(j - 1, L(sys)) + 1
    col = j % L(sys)

    println("row ", row)
    if 1 < row < L(sys)  && 1 < col < L(sys)
        onsite = bias(sys)

    else
        onsite = 0.0
    end 

    return onsite
end 




#Onsite(sys::DPT, j) = j < L(sys) ? bias_L(sys) : j < L(sys) + 1 ? bias_doubledot(sys)[1] : j < L(sys) + 2 ? bias_doubledot(sys)[2] : bias_R(sys)

function Onsite(sys::DPT, j ::Int) :: Float64


        #lower
    if j == dd_lower(sys)
        onsite = bias_doubledot(sys)[1]
        # left_offset from int term
        onsite -= U(sys) *  couple_range(sys)

    #upper
    #elseif j == L(sys) + 2
    elseif j == dd_lower(sys) + 1
        onsite = bias_doubledot(sys)[2]

    elseif j <= L_end(sys)
        onsite = bias_L(sys)

        # left_offset from int
        if j >= L_contact(sys)
            onsite -= 1/2 * U(sys)
        end 

    else
        onsite = bias_R(sys)

        # left_offset from int
        if j <= R_contact(sys)
            onsite -= 1/2 * U(sys) 
        end 

    end

    return onsite

end 

# Onsite interactions for mixed,
function Onsite(sys::DPT_mixed, j::Int)

    mix_onsite(sys, j_mix) = energies(sys)[j_mix] + (LR(sys)[j_mix] > 0 ? bias_L(sys) : bias_R(sys))
    # 'true' reservoir L, regardless of contact region

    # lower
    # we do these two first in case we move the dd sites around
    if j == dd_lower(sys) 
        onsite = bias_doubledot(sys)[1]
        # left_offset from int term
        onsite -= U(sys) *  couple_range(sys)

    #upper
    elseif j == dd_lower(sys) + 1
        onsite = bias_doubledot(sys)[2]

    # left no contact
    elseif j < L_contact(sys)

        #adjust
        j_mix = j - L_begin(sys) + 1
        onsite = mix_onsite(sys, j_mix)

    # left contact region
    elseif j <= L_end(sys)

        if !includeU(sys)
            onsite = bias_L(sys) - 1/2 * U(sys)

        else
            j_mix = j - L_begin(sys) + 1
            onsite = mix_onsite(sys, j_mix)
        end 
    
    # right contact region
    elseif j <= R_contact(sys)

        if !includeU(sys)
            onsite = bias_R(sys) - 1/2 * U(sys)

        else
            j_mix = j - L_begin(sys) + 1 - (R_begin(sys) - L_end(sys)  - 1)
            onsite = mix_onsite(sys, j_mix)
        end 

    # right no contact
    else
        
        if !includeU(sys)
            j_mix = j - L_begin(sys) + 1 - (R_begin(sys) - L_end(sys) - 1) - 2 * couple_range(sys)
        else
            j_mix = j - L_begin(sys) + 1 - (R_begin(sys) - L_end(sys) - 1) 
        end 
        
        onsite = mix_onsite(sys, j_mix)

    end 

    return onsite

end 

Onsite(sys::DPT_graph, j::Int) = Onsite(sys.dpt, j)


function Onsite(sys::DPT_avg, j::Int) 

    if j == dd_lower(sys)
        return bias_L(sys) - 1/2 * U(sys), bias_doubledot(sys)[1] - U(sys) * couple_range(sys)

    elseif j == dd_lower(sys) + 1
        return bias_R(sys) - 1/2 * U(sys), bias_doubledot(sys)[2] 

    else
        
        onsite = Onsite(sys.dpt, j)
        return onsite, 0

    end 

end

function Onsite(sys::LSR_SIAM, j::Int)


    if j <= L(sys) 
        onsite = bias_L(sys)

    elseif j == L(sys) + 1
        onsite = bias_onsite(sys)

    else
        onsite = bias_R(sys)

    end 

    return onsite

end 


Onsite(sys::reservoir_spatial, j::Int; left_offset=0.0) = sys.bias


function Onsite(sys::SD_array, j::Int)

    source = get_systotal(sys.source)
    array = get_systotal(sys.array)

    if j <= source
        return Onsite(sys.source, j)

    elseif j <= source + array
        return Onsite(sys.array, j; left_offset=source)

    else
        return Onsite(sys.drain, j; left_offset= source + array)

    end 

end 


onsiteoperators(sys::systems) = systype(sys) == "Fermion" ? ["N"] : ["Ntot"]
onsiteoperators(sys::DPT_avg) = ["Nup", "Ndn"]

function add_onsite!(sys::systems, res::OpSum)
    
    @info "Adding onsite"


    for j=1 :get_systotal(sys)
        # E-E and N-

        onsite..., = Onsite(sys, j)

        for (i, op) in enumerate(onsiteoperators(sys))

            res += onsite[i], op, sitemap(sys, j)

        end 


    end 

    return res

end 
