"""we completely decoupled the code logic of static simulations, it is required that one provides an initial state"""
function run_static_simulation(sys::systems, simulation::Static, ψ::MPS; message = "Static")

    println(message)
    h = gen_hamiltonian(sys)

    #@show h 
    saveham(message, h)

    H = MPO(h, siteinds(ψ))

    solve(H, ψ, simulation)
end 

function run_static_simulation(sys::QE_G_SIAM, simulation::Static, ψ)

    h = gen_hamiltonian(sys)

    #@show h 

    saveham(message, h)

    H = ttn(h, siteinds(ψ))
    solve(H, ψ, simulation)

end 


"""we completely decoupled the code logic of dynamic simulations, it is required that one provides an initial state"""

function run_dynamic_simulation(sys::systems, simulation::Dynamic, ψ::MPS; message = "Dynamic")

    println(message)
    h = gen_hamiltonian(sys)
    H = MPO(h, siteinds(ψ))

    #@show h

    saveham(message, h)

    time_evolve(H, ψ, simulation)

end 





