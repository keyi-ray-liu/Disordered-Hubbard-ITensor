function Onsite(sys, j)

end 

function Onsite(sys::QE_two, j) :: Float64

    systotal = get_systotal(sys)
    QEen = QEen(sys)

    if j == 1 || j == systotal
        onsite = 0.0

    elseif j == 2 || j == systotal - 1
        onsite = QEen

    else
        # adjust position of j
        onsite = Onsite(sys.chain_only, j - 2)
    end

    return den
end 

function Onsite(sys::Chain_only, j) :: Float64

    _, λ_ne, _, _, range, CN, ζ = CoulombParameters(sys)

    λ_ne *= CN
    systotal = get_systotal(sys)

    onsite = - λ_ne * ( sum([ 1/(dis(j, k, sys) + ζ) for k in max(1, j - range) : min(systotal, j + range)])   - 1/ζ)


    return onsite

end 