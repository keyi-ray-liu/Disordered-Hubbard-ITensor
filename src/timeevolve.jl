"""worker function perfomring time evolution"""
function time_evolve(ψ, sites, paras, start, fin, τ)

    QE = paras["QE"]
    L = paras["L"]
    cutoff = paras["TEcutoff"]
    maxdim = paras["TEdim"]
    method = paras["TEmethod"]
    disorder = paras["disorder"]
    

    if QE == 0
        throw(ArgumentError("Dynamics have to include quantum emitters"))

    elseif QE > 2
        throw(ArgumentError("More than two QE, not allowed"))
    end 


    prefix = getworkdir()
    disx, disy = setdisorder(disorder, L)

    #currently only support one case
    disx = disx[1, :]
    disy = disy[1, :]

    if method == "TEBD"


        gates = init_ham(paras, L, disx, disy, sites, if_gate=true)
        H = init_ham(paras, L, disx, disy, sites )
        # reverse gates
        append!(gates, reverse(gates))


        for dt in start:τ:fin
            #Sz = expect(psi, "Sz"; sites=c)

            # testing the occupation on the whole chain plus emitter
            println("$method time : $dt")
        
            ψ = apply(gates, ψ; cutoff=cutoff, maxdim = maxdim)
            
            normalize!(ψ)
            var = variance(H, ψ)

            println("variance of the state: $var")

            wf = h5open( prefix * "t$method" * string(dt) * ".h5", "w")
            write(wf, "psi", ψ)
            close(wf)
            
        end

    ################
    # TDVP starts
    elseif method == "TDVP"

        # same as DMRG, we set up H first
        H = init_ham( paras, L, disx, disy, sites)
        #forward evolution

        for dt in start:τ:fin

            println("$method time : $dt")
            ψ1 = tdvp(H, ψ, -1.0im * dt; nsweeps=1, cutoff, nsite=1)

            println( "inner", abs(inner(ψ1, ψ)))
            ψ = ψ1

            var = variance(H, ψ)
            println("variance of the state: $var")


            wf = h5open( prefix * "t$method" * string(dt) * ".h5", "w")
            write(wf, "psi", ψ)
            close(wf)

        end 

    end 

    return nothing

end 