
DenDenNeighbor(sys::systems, j) = []

function DenDenNeighbor(sys::QE_flat_SIAM, j)

end 

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
            res += U, ops, j, ops, k

        end

    end 

    return res

end 

