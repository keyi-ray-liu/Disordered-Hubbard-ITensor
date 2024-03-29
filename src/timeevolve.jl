"""worker function perfomring time evolution"""
function time_evolve(ψ, sites, paras, start, fin, occ_direct; kwargs...)

    #QE = paras["QE"]
    L = paras["L"]
    cutoff = paras["TEcutoff"]
    maxdim = paras["TEdim"]
    method = paras["TEmethod"]
    disorder = paras["disorder"]
    τ = paras["τ"]
    type = paras["type"]
    adiabatic_time = paras["adiabatic_time"]
    #source_config = paras["source_config"]
    #drain_config = paras["drain_config"]
    prefix = getworkdir()
    variance_flag = get(kwargs, :variance_flag, false)

    # if QE == 0 && s_len + length(drain)
    #     throw(ArgumentError("Dynamics have to include quantum emitter/SD"))

    if start > fin && τ > 0
        throw(ArgumentError("reverse parameters not correct"))
    end 

    if start == τ
        gs = h5open( prefix * "gs.h5", "w")
        write(gs, "psi1", ψ)
        close(gs)
    end 

    disx, disy = setdisorder(disorder, L)

    #currently only support one case
    disx = disx[1, :]
    disy = disy[1, :]

    if method == "TEBD"


        gates = init_ham(paras, L, disx, disy, sites; kwargs..., if_gate=true)
        H = init_ham(paras, L, disx, disy, sites )
        

        for dt in start:τ:fin
            #Sz = expect(psi, "Sz"; sites=c)

            # testing the occupation on the whole chain plus emitter
            println("$method time : $dt")
        
            ψ = apply(gates, ψ; cutoff=cutoff, maxdim = maxdim)

            if occ_direct
                occ = expect(ψ, "N")

                open(prefix * "occ" * string(τ), "a") do f
                    writedlm( f, transpose(occ))
                end 
            end 
            
            normalize!(ψ)

            if variance_flag
                var = variance(H, ψ)
                println("variance of the state: $var")
            end 

            wf = h5open( prefix * "t$method" * string(dt) * ".h5", "w")
            write(wf, "psi", ψ)
            close(wf)
            
        end

    ################
    # TDVP starts
    elseif method == "TDVP"

        # same as DMRG, we set up the full Hamiltonian first
        H_full = init_ham( paras, L, disx, disy, sites; kwargs...)
        #forward evolution

        for dt in start:τ:fin

            if dt < adiabatic_time

                temp_para = deepcopy(paras)
                temp_para["dp"] *= dt/adiabatic_time

                print(temp_para["dp"])
                H = init_ham( temp_para, L, disx, disy, sites; kwargs...)

            else
                H = H_full
            end 
                

            println("$method time : $dt")
            #ψ1 = tdvp(H, ψ, -1.0im * τ;  nsweeps=20, cutoff, nsite=2)
            ψ1 = tdvp(H,  -im * τ, ψ; maxdim = maxdim,  cutoff=cutoff, nsite=2, time_step= -im * τ/2, normalize=true)

            println( "inner", abs(inner(ψ1, ψ)))
            ψ = ψ1
            
            if occ_direct

                if type == "Electron"
                    oup = expect(ψ, "Nup")
                    open(prefix * "occup" * string(τ), "a") do f
                        writedlm( f, transpose(oup))
                    end 
                    
                    odn = expect(ψ, "Ndn")
                    open(prefix * "occdn" * string(τ), "a") do f
                        writedlm( f, transpose(odn))
                    end 

                elseif type == "Fermion"
                    occ = expect(ψ, "N")
                    open(prefix * "occ" * string(τ), "a") do f
                        writedlm( f, transpose(occ))
                    end 
                end
            end 

            if variance_flag
                var = variance(H, ψ)
                println("variance of the state: $var")
            end 


            wf = h5open( prefix * "t$method" * string(dt) * ".h5", "w")
            write(wf, "psi", ψ)
            close(wf)

        end 
    
    elseif occursin("eigen", method)

        staticenergy, staticwf, overlaps = load_eigen(ψ)

        if method == "eigen-occ"
            res = []
            tcd_dict = load_tcd()
        end 
        
        if method == "eigen-tcd"
            tcd_dict = load_tcd()
            num = minimum( map(maximum, collect(zip(collect(keys(tcd_dict))...))))
            res_tcd = [[] for _ in 1:num]
            res_gpi = []
        end 

        for dt in start:τ:fin

            println("$method time : $dt")
            phase = exp.( -im * dt * staticenergy)
            
            if method == "eigen-wf"
                
                teeigenwf = phase .* overlaps .* staticwf
                totalwf = add(teeigenwf...; cutoff=cutoff)

                println("overlap now: ", abs(inner( totalwf', ψ)))
                wf = h5open( prefix * "t$method" * string(dt) * ".h5", "w")
                write(wf, "psi", totalwf)
                close(wf)

            elseif method == "eigen-occ"
                
                phased_overlap = phase .* overlaps
                occ = get_eigen_occ(tcd_dict, phased_overlap)
                append!(res, [occ])

            elseif method == "eigen-tcd"

                phased_overlap = phase .* overlaps
                tcd_gs = get_eigen_tcd(tcd_dict, phased_overlap)
                gpi_gs = get_gpi(tcd_gs, paras)
                
                append!(res_gpi, gpi_gs)
                for ex in 1:num
                    append!( res_tcd[ex], [tcd_gs[ex]])
                end 

            end 

        end 

        if method == "eigen-occ"
            writedlm(prefix * "occ", res)
        end 

        if method == "eigen-tcd"
            for ex in 1:num
                writedlm( prefix * "tcd_gs" * string(ex), res_tcd[ex])
            end 

            writedlm( prefix * "gpi_gs", res_gpi)
        end 

    else
        error(ArgumentError("Unrecognized method"))

    end 

    return nothing

end 



