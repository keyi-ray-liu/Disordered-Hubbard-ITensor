"""worker function perfomring time evolution"""
function time_evolve(ψ, sites, paras, start, fin, τ)


    cutoff = paras["TEcutoff"]
    maxdim = paras["TEdim"]
    method = paras["TEmethod"]
    disorder = para["disorder"]
    

    if QE == 0
        throw(ArgumentError("Dynamics have to include quantum emitters"))

    elseif QE > 2
        throw(ArgumentError("More than two QE, not allowed"))
    end 


    prefix = getworkdir()
    disx, disy = setdisorder(disorder, L)

    if method == "TEBD"


        gates = init_ham(paras, L, disx, disy, sites, if_gate=true)
        # reverse gates
        append!(gates, reverse(gates))


        for dt in start:τ:fin
            #Sz = expect(psi, "Sz"; sites=c)

            # testing the occupation on the whole chain plus emitter
            println("$dt")
        
            ψ = apply(gates, ψ; cutoff=cutoff, maxdim = maxdim)
            normalize!(ψ)

            wf = h5open( prefix * "t" * string(dt) * ".h5", "w")
            write(wf, "psi", ψ)
            close(wf)
            
        end

    ################
    # TDVP starts
    elseif method == "TDVP"

        # same as DMRG, we set up H first
        H = init_ham( paras, L, disx, disy, sites)
        tdvp(H, ψ, t* 2im)

    end 

    return nothing

end 