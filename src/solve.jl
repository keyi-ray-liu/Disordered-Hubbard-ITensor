
""" DMRG routine that solves for the states given H and initial guess.

Returns: array of MPS's"""
function solve(H::MPO, ϕ::MPS, simulation::StaticSimulation, workflag :: String ) 

    ex, prev_state, prev_energy, prev_var, sweepcnt, sweepdim, noise, TEcutoff, krylovdim, weight = SimulationParameters(simulation)

    prev_state = convert(Vector{MPS}, prev_state)

    workdir = getworkdir(workflag)
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
            

            energy, ψ = dmrg(H, ϕ, sweeps; eigsolve_krylovdim = krylovdim, #observer=SizeObserver()
            )
            cur_ex += 1
            

        else

            energy, ψ = dmrg(H, prev_state, ϕ, sweeps; weight=weight, eigsolve_krylovdim = krylovdim, 
            #observer=SizeObserver()
            ) 

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
        wf = h5open( workdir * TEMP_tag * out * string((cur_ex - 1)) * ".h5", "w")
        write(wf, "psi", ψ)
        close(wf)

        # shift and invert block

        open( workdir * "tempex_" * out, "a") do io
            writedlm(io, energy)
        end


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

            try
                writedlm(io, expect(psi, "N"))
            catch 
                writedlm(io, expect(psi, "Ntot"))
            end 
            
        end
    end


    return prev_state
    #return prev_energy, prev_state, allenergy, allvars, prev_var

    

end 


function time_evolve(H::MPO, ψ::MPS, simulation::DynamicSimulation, workflag ; save_every=true, obs=Function[], corr_cutoff=Inf, init_obs = true, kwargs...)

    τ, start, fin, TEcutoff, TEdim, nsite= SimulationParameters(simulation)

    # get the t=0 stats

    if init_obs
        for ob in obs
            ob(;ψ=ψ, t=0, corr_cutoff = corr_cutoff,  workflag = workflag, kwargs...)
        end 

        open(getworkdir(workflag) * "times", "a" ) do io
            writedlm(io, start - τ)
        end 
    end 

    for dt in start:τ:fin

        # step(; sweep) = sweep
        # state_size(;state) = Base.format_bytes(Base.summarysize(state))
        # sys_obs = observer(
        #   "sizes" => state_size, "steps" => step
        # )


        @info "TDVP time : $dt"
        #ψ1 = tdvp(H, ψ, -1.0im * τ;  nsweeps=20, TEcutoff, nsite=2)
        @time ψ = tdvp(H,  -im * τ, ψ; updater_kwargs = (; eager=true), maxdim = TEdim,  cutoff=TEcutoff, nsite=nsite, time_step= -im * τ, normalize=true,outputlevel=1,
        # (step_observer!)=sys_obs
        )

        mem = Sys.maxrss()/2^30
        @info "Max. RSS = $mem GB"

        # println("\nResults")
        # println("=======")
        # for n in 1:length(sys_obs.steps)
        #   println("step = ", sys_obs.steps[n])
        #   println("After sweep |psi| =", sys_obs.sizes[n] )
        # end

        #println( "inner", abs(inner(ψ1, ψ)))
        #ψ = ψ1

        # we might need to calculate observables on the go
        for ob in obs
            ob(;ψ=ψ, t=dt, corr_cutoff = corr_cutoff, workflag = workflag, kwargs...)
        end 

        if save_every
            wfstr = getworkdir(workflag) * "tTDVP" * string(dt) * ".h5"
        else
            wfstr = getworkdir(workflag) * "tTDVPlaststate.h5"

            open(getworkdir(workflag) * "tTDVPlasttime", "w") do io
                writedlm(io, dt)
            end

        end 

        open(getworkdir(workflag) * "times", "a" ) do io
            writedlm(io, dt)
        end 

        h5open(wfstr, "w") do io
            write(io, "psi1", ψ)
        end 


    end 

    return ψ

end 

