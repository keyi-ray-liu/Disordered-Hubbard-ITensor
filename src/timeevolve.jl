"""worker function perfomring time evolution"""
function time_evolve(ψ; t=1.0, L=12, N=6, CN=6, QE=2, QN=false, QEen=0.4, dp=0.0075, time=5, τ=0.1, disorder=false)

    if QE == 0
        throw(ArgumentError("Dynamics have to include quantum emitters"))

    elseif QE > 2
        throw(ArgumentError("More than two QE, not allowed"))
    end 

    # get disorder (if no then zeros) 
    disx, disy = setdisorder(disorder, L)

    sites =  siteinds("Fermion", N + QE; conserve_qns=true)
    gates = ITensor[]
    
    # two-site hopping gates
    for j in 1:(L - 1)

        r = dis(j, j + 1, disx, disy)
        hop = hopping(decay, r)

        s1 = sites[j]
        s2 = siste[j + 1]
        hj =
          - t * hop * op("C", s1) * op("Cdag", s2) +
          - t * hop * op("C", s2) * op("Cdag", s1)
        Gj = exp(-im * τ / 2 * hj)
        push!(gates, Gj)
    end

    # reverse gates
    append!(gates, reverse(gates))

    # adding hopping gates 

end 