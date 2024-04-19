"""we completely decoupled the code logic of static simulations, it is required that one provides an initial state"""
function run_static_simulation(sys::systems, simulation::Static, ψ::MPS)

    h = gen_hamiltonian(sys)

    #@show h 
    H = MPO(h, siteinds(ψ))

    solve(H, ψ, simulation)
end 

function run_static_simulation(sys::QE_G_SIAM, simulation::Static, ψ)

    h = gen_hamiltonian(sys)

    @show h 
    H = ttn(h, siteinds(ψ))
    solve(H, ψ, simulation)

end 


"""we completely decoupled the code logic of dynamic simulations, it is required that one provides an initial state"""

function run_dynamic_simulation(sys::systems, simulation::Dynamic, ψ::MPS)

    h = gen_hamiltonian(sys)
    H = MPO(h, siteinds(ψ))

    @show h
    time_evolve(H, ψ, simulation)

end 


function run_Static_chain(dims, N;  kwargs...)


    coulomb = set_Coulombic(;kwargs...)
    sys = set_Chain(; dims=dims, N=N, coulomb=coulomb, kwargs...)

    simulation = set_Static(;kwargs...)
    ψ = gen_state(sys)

    run_static_simulation(sys, simulation, ψ)

end 


