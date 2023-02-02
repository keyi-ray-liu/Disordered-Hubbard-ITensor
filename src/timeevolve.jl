"""worker function perfomring time evolution"""
function time_evolve(ψ; t=1.0, L=12, N=6, CN=6, QE=2, QN=false, QEen=0.4, dp=0.0075, t=5, τ=0.1)

    if QE == 0
        throw(ArgumentError("Dynamics have to include quantum emitters"))

    elseif QE > 2
        throw(ArgumentError("More than two QE, not allowed"))
    end 

    sites =  siteinds("Fermion", N + QE; conserve_qns=true)
    gates = ITensor[]

    # adding hopping gates 

end 