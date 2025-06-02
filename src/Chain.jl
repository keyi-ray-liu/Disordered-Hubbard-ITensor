function run_GQS(L, N; τ=1.0, fin=1000.0, kwargs...)

    sys = set_GQS(; L=L, N=N, kwargs...)
    simulation = DynamicSimulation(;start=τ, fin=fin, τ=τ, kwargs...)
    ψ = gen_state(sys)

    run_dynamic_simulation(sys, simulation, ψ)

    # ψ = load_ψ("wf")
    # partial_contract(ψ, [4, 5])
end 


function run_rectangular()

    sys = Rectangular(; Lx = 3, Ly = 3, 
    N = 4, systype = "Fermion", 
    #N = [2, 3, 4, 0], systype = "Electron", 
    U = 4.0, 
    )

    static = StaticSimulation(; ex=1, sweepdim=400)
    ψ = gen_state(sys)
    ψ = run_static_simulation(sys, static, ψ, Identity())

    return ψ
    

end 





function run_chain( template::Union{typeof(Chain), typeof(Ring)}, L, N, ex; dim=64, kwargs...)

    sys = template(;  L=L, N=N, kwargs...)

    @show sys
    static = StaticSimulation(; ex=ex, sweepdim=dim, kwargs...)
    ψ = gen_state(sys)
    ψ = run_static_simulation(sys, static, ψ, Identity())

    return ψ

end 

function run_ring( L, N, ex; dim=64, kwargs...)

    sys = Chain(;  L=L, N=N, kwargs...)

    @show sys
    static = StaticSimulation(; ex=ex, sweepdim=dim, kwargs...)
    ψ = gen_state(sys)
    ψ = run_static_simulation(sys, static, ψ, Identity())

    return ψ

end 


function run_biased_chain(fullsize, L, N, ex; dim=64, sweepcnt=40, kwargs...)

    sys = set_biased_chain(; biaswindow = [1, L], L=fullsize, N=N, kwargs...)

    @show sys
    static = StaticSimulation(; ex=ex, sweepdim=dim, sweepcnt=sweepcnt, kwargs...)
    ψ = gen_state(sys; sites=get(kwargs, :sites, nothing))
    ψ = run_static_simulation(sys, static, ψ, Identity())

    return ψ

end 

function biased_quench(L, N; biaswindow=[1, 10], bias=500, dim=64, τ=0.5,  tswitch=50, fin=100, kwargs...)

    sys = Chain(;  L=L, N=N, kwargs...)
    obs = [dyna_EE, dyna_occ, dyna_corr]

    @show sys
    @info "Solve GS"
    static = StaticSimulation(; ex=1, sweepdim=dim, sweepcnt=10, kwargs...)
    ψ = gen_state(sys)
    ψ = run_static_simulation(sys, static, ψ, Identity())[1]

    @info "Quench"
    biased_sys = set_biased_chain(; L=L, N=N, biaswindow=biaswindow, bias=bias)
    @show biased_sys
    Quench =  Dynamic(; TEdim=dim, τ=τ, start= τ, fin=tswitch, kwargs...)
    ψ = run_dynamic_simulation(biased_sys, Quench, ψ; message="Quench", save_every=false, obs=obs)
    
    @info "Relaxation"
    Relaxation = DynamicSimulation(; TEdim=dim, τ=τ, start= tswitch + τ, fin=fin, kwargs...)
    _ = run_dynamic_simulation(sys, Relaxation, ψ; message="Relaxation", save_every=false, obs=obs)



    return nothing

end 



function GQS_wrapper()

    chain_in = load_JSON(pwd() * "/GQS.json")

    L = get(chain_in, "L", 12)
    N = get(chain_in, "N", 6)
    sweepcnt = get(chain_in, "sweepcnt", 20)
    init = get(chain_in, "init", 1)

    λ_ee = get(chain_in, "int_ee", 3.0)
    λ_ne = get(chain_in, "int_ne", 3.0)
    

    run_GQS(L, N; sweepcnt=sweepcnt, init=init, λ_ee =λ_ee, λ_ne = λ_ne)
    

end 


function chain_wrapper()

    chain_in = load_JSON(pwd() * "/chain.json")

    L = get(chain_in, "L", 12)
    N = get(chain_in, "N", 6)
    sweepcnt = get(chain_in, "sweepcnt", 20)
    ex = get(chain_in, "ex", 10)
    dim = get(chain_in, "dim", 64)
    

    run_chain( Chain, L, N, ex; sweepcnt=sweepcnt, dim=dim)

    return nothing
end 



# function rectangular_wrapper()


#     rect_in = load_JSON(pwd() * "/rectagular.json")

#     Lx = get(rect_in, "Lx", 3)
#     Ly = get(rect_in, "Ly", 3)
#     N = get(rect_in, "N", 0)
#     sweepcnt = get(rect_in, "sweepcnt", 20)
#     ex = get(rect_in, "ex", 1)
#     dim = get(rect_in, "dim", 64)
    

#     run_chain( Chain, L, N, ex; sweepcnt=sweepcnt, dim=dim)

#     return nothing
# end 



function quench_wrapper()

    chain_in = load_JSON(pwd() * "/chainquench.json")

    L = get(chain_in, "L", 12)
    N = get(chain_in, "N", 6)
    
    τ = get(chain_in, "timestep", 0.5)
    tswitch = get(chain_in, "tswitch", 50)
    fin = get(chain_in, "fin", 100)
    biaswindow = get(chain_in, "biaswindow", [1, 10])
    dim = get(chain_in, "dim", 64)
    bias = get(chain_in, "bias", 500)

    biased_quench(L, N; biaswindow = biaswindow,dim=dim, τ=τ, tswitch=tswitch, fin=fin, bias=bias)

    return nothing
end 


function ring_wrapper()

    ring_in = load_JSON(pwd() * "/ring.json")

    L = get(ring_in, "L", 12)
    N = get(ring_in, "N", 6)
    sweepcnt = get(ring_in, "sweepcnt", 20)
    systype = get(ring_in, "systype", "Electron")
    ex = get(ring_in, "ex", 1)
    dim = get(ring_in, "dim", 128)
    U = get(ring_in, "U", 4.0)
    

    run_chain(Ring(), L, N, ex; sweepcnt=sweepcnt, dim=dim, U=U, systype=systype)

    return nothing
end 