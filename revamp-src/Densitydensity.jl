


function DenDenNeighbor(sys::QE_flat_SIAM, j)

end 

function DenDenNeighbor(sys::QE_two, j)

    systotal = get_systotal(sys)

    if j <3 || j > systotal - 3
        den = []

    else
        # shift the position of each j
        den = DenDenNeighbor(sys.chain_only, j - 2)

        # shift back position of each k
        den = [ [U, k + 2] for (U, k) in den]
    end

    return den
end 

function DenDenNeighbor(sys::Chain_only, j)

    λ_ee, _, exch, _, range, _, ζ = CoulombParameters(sys)
    ifexch(j, k, sys) = ( 1 - (dis(j, k, sys) == 1) * exch )
    
    return [ [λ_ee * ifexch(j, k, sys) / ( dis(j, k, sys) + ζ), k] for k in max(1, j - range) : j - 1]
    
end 