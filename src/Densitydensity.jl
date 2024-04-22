

""" A general principle is that we always use native indexing in the subsystem"""
DenDenNeighbor(sys::systems, j) = []

function DenDenNeighbor(sys::QE_flat_SIAM, j)


    _, chain_begin, chain_end = get_sys_loc(sys, j)
    λ_ee, _, exch, _, range, _, ζ = CoulombParameters(sys)

    ifexch(j, k, sys) = ( 1 - (dis(j, k, sys) == 1) * exch )

    if chain_begin <= j <= chain_end

        return [ [λ_ee * ifexch(j, k, sys) / ( dis(j, k, sys) + ζ), k] for k in max(chain_begin, j - range) : j - 1]


    # center site interacting, since the sys is completely symmetric, the sites are the only thing thats being changed
    elseif j == left(sys) + 1

        denden = []

        for k in 1:min(range, siteseach(sys))

            # we calculate interact ONCE, as if we are at site 0
            strength = center_ee(sys) * ifexch(0, k, sys) / (dis(0, k, sys) + ζ)
            #left
            for l in 1:legleft(sys)
                append!(denden, [[strength, j - k - (l -1) * (siteseach(sys) + QESITES)]])
            end 

            #right
            for r in 1:legright(sys)
                append!(denden, [[strength, j + k + (r -1) * (siteseach(sys) + QESITES)]])
            end 

        end 

        @show denden
        return denden
    else

        return []
    end 

end 

DenDenNeighbor(sys::QE_G_SIAM, j) = DenDenNeighbor(sys.system, j)

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

function DenDenNeighbor(sys::QE_parallel, j)

    uppertotal = get_uppertotal(sys)
    systotal = get_systotal(sys)

    if j <= uppertotal
        den =  DenDenNeighbor(sys.upper, j)

    elseif j < systotal
        den =  [ [U, k + uppertotal] for (U, k) in DenDenNeighbor(sys.lower, j - uppertotal)]

    # the contact site, physically sits in the middle of two chains. We move to the native chain coordinates
    else
        upperchain = get_upperchain(sys)
        lowerchain = get_lowerchain(sys)

        uppercenter = div(upperchain, 2) + 1
        lowercenter = div(lowerchain, 2) + 1
        
        upperoffset = QESITES
        loweroffset = uppertotal + QESITES

        den = vcat(
            [ [ center_ee(sys) / (dis(i, uppercenter, sys) + center_dis(sys)), i + upperoffset] for i in  max( 1, uppercenter - center_range(sys) ) : min( upperchain, uppercenter + center_range(sys))],
            [ [ center_ee(sys) / (dis(i, lowercenter, sys) + center_dis(sys)), i + loweroffset] for i in  max( 1, lowercenter - center_range(sys) ) : min( lowerchain, lowercenter + center_range(sys))]
        )
    end 

    return den
end 

function DenDenNeighbor(sys::Chain_only, j::Int)

    λ_ee, _, exch, _, range, _, ζ = CoulombParameters(sys)
    ifexch(j, k, sys) = ( 1 - (dis(j, k, sys) == 1) * exch )
    
    return [ [λ_ee * ifexch(j, k, sys) / ( dis(j, k, sys) + ζ), k] for k in max(1, j - range) : j - 1]
    
end 

function DenDenNeighbor(sys::DPT, j::Int)

    # if j == L(sys) + 1
    if j == L(sys) + R(sys) + 1

        #den =  vcat( [ [U(sys), j - k] for k in 1:couple_range(sys)], [[U(sys), j + 1 + k] for k in 1:couple_range(sys)])

        den =  vcat( [ [U(sys), L(sys) + 1 - k] for k in 1:couple_range(sys)], [[U(sys), L(sys) + k] for k in 1:couple_range(sys)])
    else
        den = []
    end 

    return den

end 


function add_DensityDensity!(sys::systems, res::OpSum)
    

    println("Adding EE and NE")
    sys_type = type(sys)

    systotal = get_systotal(sys)

    if sys_type == "Fermion"
        ops = "N"
  
    elseif sys_type == "Electron"
        ops = "Ntot"
  
    end 

    for j=1 :systotal
        # E-E and N-E

        for (U, k) in DenDenNeighbor(sys, j)
            
            k = trunc(Int, k)
            # delta function setting up the exchange between nearest neighbor

            # because k < j, we check if j is in the nn of k
            res += U, ops, sitemap(sys, j), ops, sitemap(sys, k)

        end

    end 

    return res

end 

