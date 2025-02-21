"""we completely decoupled the code logic of static SimulationParameters, it is required that one provides an initial state
Returns: return of solve function. Array of MPS's
"""
function run_static_simulation(sys::Systems, simulation::StaticSimulation, ψ::MPS; process :: StateModifier = Identity(), message = "Static")

    @info message
    @show sys
    h = gen_hamiltonian(sys)

    #@show h 
    saveham(message, h)

    H = MPO(h, siteinds(ψ))

    state = solve(H, ψ, simulation)
    state = [modifystate(st, process, sys) for st in state]
    return state
end 


"""we completely decoupled the code logic of dynamic SimulationParameters, it is required that one provides an initial state"""

function run_dynamic_simulation(sys::Systems, simulation::DynamicSimulation, ψ::MPS; message = "Dynamic", save_every=true, obs=Function[], init_obs=true,  kwargs...)

    @info message
    h = gen_hamiltonian(sys)
    H = MPO(h, siteinds(ψ))

    #@show h

    saveham(message, h)

    ψ = time_evolve(H, ψ, simulation; save_every=save_every, obs=obs, sys=sys, init_obs=init_obs, kwargs...)

    return ψ

end 

"""Wrapper function to automatically load last checkpoint, and run simulations according to the given parameters and system configurations"""
function run_gs_dyna(timecontrol :: OneStage, init::Union{Nothing, Systems}, sys::Systems, obs; random=true, process :: StateModifier = Identity(), kwargs...)

    τ = timecontrol.τ
    fin = timecontrol.fin

    last_time, last_state = prev_res()
    @show last_time
    start = last_time > -Inf ? last_time + τ : τ

    ψ = gen_state(init; random=random, kwargs...)
    # we solve for the GS of the whole system at zero bias, we also bias the array so that nothing is occupied there

    if isnothing(init)
        ψ0 = ψ
    else
        Static = StaticSimulation(; output=EQINIT_STR, sweepdim=get(kwargs, :TEdim, 256) * 2, sweepcnt=get(kwargs, :sweepcnt, 50), ex=1, kwargs...)

        # GS calculation
        ψ0 :: MPS =  last_time > -Inf ? last_state : check_ψ(EQINIT_STR) ? load_ψ(EQINIT_STR) : run_static_simulation(init, Static, ψ; message = "Init", process = process)[1]

    end 

    Dynamic = DynamicSimulation(;τ=τ, start=start, fin=fin, kwargs...)
    _ = run_dynamic_simulation(sys, Dynamic, ψ0; save_every=false, obs=obs, kwargs...)
end 

"""Wrapper function to automatically load last checkpoint, and run a two-stage dynamical simulation, with specified t1 and t2"""
function run_gs_dyna(timecontrol::TwoStage, init::Union{Nothing, Systems}, sys::Systems, obs; random=true, process :: StateModifier = Identity(), kwargs...)

    τ1 = timecontrol.τ1
    fin1 = timecontrol.fin1
    τ2 = timecontrol.τ2
    fin2 = timecontrol.fin2

    last_time, last_state = prev_res()
    @show last_time


    ψ = gen_state(init; random=random, kwargs...)
    # we solve for the GS of the whole system at zero bias, we also bias the array so that nothing is occupied there

    if isnothing(init)
        ψ0 = ψ
    else
        Static = StaticSimulation(; output=EQINIT_STR, sweepdim=get(kwargs, :TEdim, 256) * 2, sweepcnt=get(kwargs, :sweepcnt, 50), ex=1, kwargs...)

        # GS calculation
        ψ0 =  last_time > -Inf ? last_state : check_ψ(EQINIT_STR) ? load_ψ(EQINIT_STR) : run_static_simulation(init, Static, ψ; message = "Init", process = process)[1]
    end 

    if last_time < fin1
        @info "Stage 1"
        s1start = max(0, last_time) + τ1
        Stage1 = DynamicSimulation(;τ=τ1, start=s1start, fin=fin1, kwargs...)
        ψ0 = run_dynamic_simulation(sys, Stage1, ψ0; save_every=false, obs=obs, kwargs...)
    end 

    @info "Stage 2"
    s2start = max(last_time, fin1) + τ2
    Stage2 = DynamicSimulation(;τ=τ2, start=s2start, fin=fin2, kwargs...)
    ψ0 = run_dynamic_simulation(sys, Stage2, ψ0; save_every=false, obs=obs, kwargs...)


end 




# function run_static_simulation(sys::Union{QE_G_SIAM, DPT_graph}, simulation::StaticSimulation, ψ; message = "Static")

#     @info "GRAPH " * message
#     h = gen_hamiltonian(sys)

#     #@show h 

#     saveham(message, h)

#     H = ttn(h, siteinds(ψ))
#     return solve(H, ψ, simulation)

# end 


# function run_dynamic_simulation(sys::Union{QE_G_SIAM, DPT_graph}, simulation::DynamicSimulation, ψ; message = "Dynamic", save_every=true, obs=Function[])

#     @info "GRAPH " * message
#     h = gen_hamiltonian(sys)

#     #@show h 
#     saveham(message, h)

#     H = ttn(h, siteinds(ψ))
#     ψ = time_evolve(H, ψ, simulation; save_every=save_every, obs=obs, sys=sys)

#     return ψ

# end 



