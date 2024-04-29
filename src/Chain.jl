function run_GQS(L, N; τ=1.0, fin=1000.0, kwargs...)

    sys = set_GQS(; L=L, N=N, kwargs...)
    simulation = set_Dynamic(;start=τ, fin=fin, τ=τ, kwargs...)
    ψ = gen_state(sys)

    run_dynamic_simulation(sys, simulation, ψ)

    # ψ = load_ψ("wf")
    # partial_contract(ψ, [4, 5])
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