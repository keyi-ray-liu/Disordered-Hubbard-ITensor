"""we completely decoupled the code logic of static SimulationParameters, it is required that one provides an initial state
Returns: return of solve function. Array of MPS's
"""
function run_static_simulation(sys::Systems, simulation::Static, ψ::MPS; message = "Static")

    @info message
    h = gen_hamiltonian(sys)

    #@show h 
    saveham(message, h)

    H = MPO(h, siteinds(ψ))

    return solve(H, ψ, simulation)
end 


"""we completely decoupled the code logic of dynamic SimulationParameters, it is required that one provides an initial state"""

function run_dynamic_simulation(sys::Systems, simulation::Dynamic, ψ::MPS; message = "Dynamic", save_every=true, obs=Function[], init_obs=true)

    @info message
    h = gen_hamiltonian(sys)
    H = MPO(h, siteinds(ψ))

    #@show h

    saveham(message, h)

    ψ = time_evolve(H, ψ, simulation; save_every=save_every, obs=obs, sys=sys, init_obs=init_obs)

    return ψ

end 

"""Wrapper function to automatically load last checkpoint, and run simulations according to the given parameters and system configurations"""
function run_gs_dyna(τ, fin, init::Union{Nothing, Systems}, sys::Systems, obs; random=true, kwargs...)

    last_time, last_state = prev_res()
    @show last_time
    start = last_time > -Inf ? last_time + τ : τ

    ψ = gen_state(sys; random=random)
    # we solve for the GS of the whole system at zero bias, we also bias the array so that nothing is occupied there

    if isnothing(init)
        ψ0 = ψ
    else
        Static = set_Static(; output=EQINIT_STR, sweepdim=get(kwargs, :TEdim, 64), sweepcnt=get(kwargs, :sweepcnt, 80), ex=1, kwargs...)

        # GS calculation
        ψ0 =  last_time > -Inf ? last_state : check_ψ(EQINIT_STR) ? load_ψ(EQINIT_STR) : run_static_simulation(init, Static, ψ; message = "Init")[1]
    end 

    Dynamic = set_Dynamic(;τ=τ, start=start, fin=fin, kwargs...)
    _ = run_dynamic_simulation(sys, Dynamic, ψ0; save_every=false, obs=obs)
end 



# function run_static_simulation(sys::Union{QE_G_SIAM, DPT_graph}, simulation::Static, ψ; message = "Static")

#     @info "GRAPH " * message
#     h = gen_hamiltonian(sys)

#     #@show h 

#     saveham(message, h)

#     H = ttn(h, siteinds(ψ))
#     return solve(H, ψ, simulation)

# end 


# function run_dynamic_simulation(sys::Union{QE_G_SIAM, DPT_graph}, simulation::Dynamic, ψ; message = "Dynamic", save_every=true, obs=Function[])

#     @info "GRAPH " * message
#     h = gen_hamiltonian(sys)

#     #@show h 
#     saveham(message, h)

#     H = ttn(h, siteinds(ψ))
#     ψ = time_evolve(H, ψ, simulation; save_every=save_every, obs=obs, sys=sys)

#     return ψ

# end 



