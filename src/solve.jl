
""" DMRG routine that solves for the states given H and initial guess.

Returns: array of MPS's"""
function solve(H::MPO, ϕ::MPS, simulation::Static) 

    ex, prev_state, prev_energy, prev_var, sweepcnt, sweepdim, noise, TEcutoff, krylovdim, weight = SimulationParameters(simulation)

    prev_state = convert(Vector{MPS}, prev_state)

    workdir = getworkdir()
    allvars = copy(prev_var)
    allenergy = copy(prev_energy)
    cur_ex = 1 + length(prev_state)
    cnt = 1
    out = output(simulation)
  
    # changed to while loop for possible retracting operation
    # global cnt to stop excessive looping
    while cur_ex <= ex && cnt <= GLOBAL_CNT

        cnt += 1
        println("len energy: ", length(prev_energy), "len state: ", length(prev_state), "cur: ", cur_ex)
        println("prev_energy so far", allenergy)
        println("variances so far", allvars)

        
        # DMRG block
        # DMRG parameters
        sweeps = Sweeps(sweepcnt)
        setmaxdim!(sweeps, sweepdim)

        if noise
            setnoise!(sweeps, 1E-5)
        end 

        setcutoff!(sweeps, TEcutoff)

        # DMRG method
        # Run the DMRG algorithm, returning energy
        # (dominant eigenvalue) and optimized MPS

        if cur_ex == 1
            

            energy, ψ = dmrg(H, ϕ, sweeps; eigsolve_krylovdim = krylovdim)
            cur_ex += 1
            

        else

            energy, ψ = dmrg(H, prev_state, ϕ, sweeps; weight=weight, eigsolve_krylovdim = krylovdim) 

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
        wf = h5open( workdir * "temp_" * out * string((cur_ex - 1)) * ".h5", "w")
        write(wf, "psi", ψ)
        close(wf)

        # shift and invert block


        println("As of end of search, systype of prev_state", typeof(prev_state))
        var = variance(H, ψ)

        append!(allenergy, energy)
        append!(prev_energy, energy)
        append!(prev_state, [ψ])
        append!(allvars, var)
        append!(prev_var, var)

    end 
    
    # save results
    

    print(prev_energy)

    writedlm( workdir * out * "ex" , prev_energy)
    writedlm( workdir * out * "var", prev_var)
    writedlm( workdir * out * "allvar", allvars)
    writedlm( workdir * out * "allenergy", allenergy)
  
    # write wf
    
    h5open(workdir * out * ".h5", "w") do io
        for (i, psi) in enumerate(prev_state)
            write(io, "psi" * string(i), psi)
        end
    end 

    open(workdir * "staticocc", "w") do io
        for psi in prev_state
            writedlm(io, expect(psi, "N"))
        end
    end


    return prev_state
    #return prev_energy, prev_state, allenergy, allvars, prev_var

    

end 


function time_evolve(H::MPO, ψ::MPS, simulation::Dynamic; save_every=true, obs=Function[], init_obs = true, kwargs...)

    # get the t=0 stats

    if init_obs
        for ob in obs
            ob(;ψ=ψ, t=0, kwargs...)
        end 
    end 

    τ, start, fin, TEcutoff, TEdim, nsite= SimulationParameters(simulation)

    for dt in start:τ:fin

        @info "TDVP time : $dt"
        #ψ1 = tdvp(H, ψ, -1.0im * τ;  nsweeps=20, TEcutoff, nsite=2)
        ψ1 = tdvp(H,  -im * τ, ψ; maxdim = TEdim,  cutoff=TEcutoff, nsite=nsite, time_step= -im * τ/2, normalize=true)

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

        end 

        open(getworkdir() * "times", "a" ) do io
            writedlm(io, dt)
        end 

        write(wf, "psi1", ψ)
        close(wf)

    end 

    return ψ

end 

