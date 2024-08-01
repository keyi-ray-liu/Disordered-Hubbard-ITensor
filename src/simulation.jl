"""we completely decoupled the code logic of static simulations, it is required that one provides an initial state"""
function run_static_simulation(sys::systems, simulation::Static, ψ::MPS; message = "Static")

    @info message
    h = gen_hamiltonian(sys)

    #@show h 
    saveham(message, h)

    H = MPO(h, siteinds(ψ))

    return solve(H, ψ, simulation)
end 


"""we completely decoupled the code logic of dynamic simulations, it is required that one provides an initial state"""

function run_dynamic_simulation(sys::systems, simulation::Dynamic, ψ::MPS; message = "Dynamic", save_every=true, obs=Function[])

    @info message
    h = gen_hamiltonian(sys)
    H = MPO(h, siteinds(ψ))

    #@show h

    saveham(message, h)

    ψ = time_evolve(H, ψ, simulation; save_every=save_every, obs=obs, sys=sys)

    return ψ

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



