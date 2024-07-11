

function solve(H::ITensorNetworks.TreeTensorNetwork{Any}, ϕ::ITensorNetworks.TreeTensorNetwork{Any}, simulation::Static) 

    ex, prev_state, prev_energy, prev_var, sweepcnt, sweepdim, noise, TEcutoff, krylovdim, weight = SimulationParameters(simulation)

    workdir = getworkdir()
    allvars = copy(prev_var)
    allenergy = copy(prev_energy)
    cur_ex = 1 + length(prev_state)
    cnt = 1
  
    # changed to while loop for possible retracting operation
    # global cnt to stop excessive looping
    while cur_ex <= ex && cnt <= GLOBAL_CNT

        cnt += 1
        println("len energy: ", length(prev_energy), "len state: ", length(prev_state), "cur: ", cur_ex)
        println("prev_energy so far", allenergy)
        println("variances so far", allvars)


        # DMRG method
        # Run the DMRG algorithm, returning energy
        # (dominant eigenvalue) and optimized MPS

        if cur_ex == 1


            #energy, ψ = dmrg(H, ϕ, sweeps; eigsolve_krylovdim = krylovdim)
            energy, ψ= dmrg(H, ϕ; nsweeps = sweepcnt, maxdim = sweepdim, cutoff = TEcutoff, nsites = 2, updater_kwargs=(; krylovdim=krylovdim), outputlevel=2)

            
            cur_ex += 1
            

        else
            energy, ψ = dmrg(H, prev_state, ϕ, sweeps; weight=weight, eigsolve_krylovdim = krylovdim, outputlevel=2) 

            # check if cur energy is lower than previously achieved energy, if so, return to the point with lower energy (if not, start with current state as GS)
            if abs(energy - prev_energy[end]) > TOL && energy < prev_energy[end]

                cur = 1

                while cur <= length(prev_energy)

                    if abs(energy - prev_energy[cur]) > TOL && prev_energy[cur] < energy
                        cur += 1
                    else 
                        break
                    end 

                end 

                # reset the current array of prev_state and prev_energy, reset ex count
                prev_energy = prev_energy[begin:cur-1]
                prev_state = prev_state[begin:cur-1]
                prev_var = prev_var[begin:cur-1]
                cur_ex = cur
            
            # else continue evaluation
            else
                cur_ex += 1

            end

            # if this is new lowest state, we reset our calculation using ψ as GS
            if cur_ex == 1

                ϕ = ψ
                cur_ex = 2

            end 
            
        end

            
        #println(expect(ψ, "N"))
        
        # save temp results
        # wf = h5open( workdir * "temp_cur_ex" * string((cur_ex - 1)) * ".h5", "w")
        # write(wf, "psi", ψ)
        # close(wf)

        # shift and invert block


        println("As of end of search, type of prev_state", typeof(prev_state))

        append!(allenergy, energy)
        append!(prev_energy, energy)
        append!(prev_state, [ψ])
        # append!(allvars, var)
        # append!(prev_var, var)

    end 
    
    # save results
    out = output(simulation)

    print(prev_energy)

    writedlm( workdir * out * "ex" , prev_energy)
    writedlm( workdir * out * "var", prev_var)
    writedlm( workdir * out * "allvar", allvars)
    writedlm( workdir * out * "allenergy", allenergy)
  
    # write wf
  
    # wf = h5open( workdir * out * ".h5", "w")
    
    # for (i, psi) in enumerate(prev_state)
    #   write(wf, "psi" * string(i), psi)
    # end
  
    # close(wf)

    #return prev_energy, prev_state, allenergy, allvars, prev_var
    return prev_state
    

end 


function time_evolve(H::ITensorNetworks.TreeTensorNetwork{Any}, ψ::ITensorNetworks.TreeTensorNetwork{Any}, simulation::Dynamic; save_every=true, obs=Function[], kwargs...)

    τ, start, fin, TEcutoff, TEdim, nsites = SimulationParameters(simulation)

    for dt in start:τ:fin


        @info "TDVP time : $dt"
        #ψ1 = tdvp(H, ψ, -1.0im * τ;  nsweeps=20, TEcutoff, nsite=2)
        ψ1 = ITensorNetworks.tdvp(H,  -im * τ, ψ; maxdim = TEdim,  cutoff=TEcutoff, nsites=nsites, time_step= -im * τ/2, normalize=true)

        #println( "inner", abs(inner(ψ1, ψ)))
        ψ = ψ1

        # we might need to calculate observables on the go
        for ob in obs
            ob(;ψ=ψ, t=dt, kwargs...)
        end 

        if save_every
            wf = h5open( getworkdir() * "tTDVP" * string(dt) * ".h5", "w")
        else
            wf = h5open( getworkdir() * "tTDVPlaststate.h5", "w")

            open(getworkdir() * "tTDVPlasttime", "w") do io
                writedlm(io, dt)
            end 

            open(getworkdir() * "times", "a" ) do io
                writedlm(io, dt)
            end 

        end 

        write(wf, "psi1", ψ)
        close(wf)

    end 

    return ψ

end 