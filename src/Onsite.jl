Onsite(sys::systems, j) = 0.0


function Onsite(sys::QE_two, j) :: Float64

    systotal = get_systotal(sys)

    if j == 2 || j == systotal
        onsite = 0.0

    elseif j == 1 || j == systotal - 1
        onsite = QEen(sys)

    else
        # adjust position of j
        onsite = Onsite(sys.chain_only, j - 2)
    end

    return onsite
end 

function Onsite(sys::QE_parallel, j) 

    uppertotal = get_uppertotal(sys)

    if j <= uppertotal
        onsite = Onsite(sys.upper, j)

    elseif j < get_systotal(sys)
        onsite = Onsite(sys.lower, j - uppertotal)

    else
        upperchain = get_upperchain(sys)
        lowerchain = get_lowerchain(sys)

        uppercenter = div(upperchain, 2) + 0.5
        lowercenter = div(lowerchain, 2) + 0.5

        onsite =  - ( sum(
            [ center_ne(sys) / (dis(i, uppercenter, sys) + center_dis(sys)) for i in  max( 1, uppercenter - center_range(sys) ) : min( upperchain, uppercenter + center_range(sys))]) + 
            sum( [ center_ne(sys) / (dis(i, lowercenter, sys) + center_dis(sys)) for i in  max( 1, lowercenter - center_range(sys) ) : min( lowerchain, lowercenter + center_range(sys))]
        ) ) 

    end 

    return onsite
    
end 

function Onsite(sys::QE_flat_SIAM, chain_begin, chain_end, j)

    _, λ_ne, _, _, range, CN, ζ = CoulombParameters(sys)

    λ_ne *= CN

    if j == left(sys) + 1

        onsite =  sum([ - center_ne(sys) * CN / (dis(0, k, sys) + ζ) for k in 1: min(range, siteseach(sys))]) * (legleft(sys) + legright(sys))


    else
        
        onsite = - λ_ne * CN *  ( sum([ 1/(dis(j, k, sys) + ζ) for k in max(chain_begin, j - range) : min(chain_end, j + range)])   - 1/ζ)
    end 

    return onsite

end 


function Onsite(sys::QE_flat_SIAM, j) :: Float64

    qe_loc, chain_begin, chain_end = get_sys_loc(sys, j)

    if j == qe_loc + 1 
        onsite = 0.0

    elseif j == qe_loc
        onsite = QEen(sys)
    else
        # adjust position of j
        onsite = Onsite(sys, chain_begin, chain_end, j)
    end

    return onsite
end 

Onsite(sys::QE_G_SIAM, j) = Onsite(sys.system, j)

function Onsite(sys::Chain_only, j::Int) :: Float64

    _, λ_ne, _, _, range, CN, ζ = CoulombParameters(sys)

    λ_ne *= CN
    systotal = get_systotal(sys)

    onsite = - λ_ne * ( sum([ 1/(dis(j, k, sys) + ζ) for k in max(1, j - range) : min(systotal, j + range)])   - 1/ζ)


    return onsite

end 

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


    if j <= L(sys)
        onsite = bias_L(sys)

        # offset from int
        if j > L(sys) - couple_range(sys)
            onsite -= 1/2 * U(sys)
        end 
    
    #lower
    # elseif j == L(sys) + 1
    elseif j == L(sys) + R(sys) + 1
        onsite = bias_doubledot(sys)[1]

        # offset from int term
        onsite -= U(sys) *  couple_range(sys)

    #upper
    #elseif j == L(sys) + 2
    elseif j == L(sys) + R(sys) + 2
        onsite = bias_doubledot(sys)[2]

    else
        onsite = bias_R(sys)

        # offset from int
        if j <= L(sys) + couple_range(sys) 
            onsite -= 1/2 * U(sys) 
        end 

    end

    return onsite

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

function add_onsite!(sys::systems, res::OpSum)
    
    println("Adding onsite")
    sys_type = type(sys)

    systotal = get_systotal(sys)

    if sys_type == "Fermion"
        ops = "N"
  
    elseif sys_type == "Electron"
        ops = "Ntot"
  
    end 

    for j=1 :systotal
        # E-E and N-
        res += Onsite(sys, j), ops, sitemap(sys, j)

    end 

    return res

end 
