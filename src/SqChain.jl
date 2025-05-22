
function run_SqChain(L, Nup, Ndn, t;  U=4.0, bias=0.0, kwargs...)


    sys = Sq_chain(; L=L, Nup=Nup, Ndn=Ndn, t=t, U=U, bias=bias, kwargs...)

    simulation = StaticSimulation(; sweepcnt=300, sweepdim = 100, kwargs...)
    ψ = gen_state(sys)

    run_static_simulation(sys, simulation, ψ, Identity())

end 



function SQ_wrapper()

    SqChain_in = load_JSON( pwd() * "/NFpara.json")

    U = get(SqChain_in, "U", 4.0)
    L = get(SqChain_in, "L", 3)
    Nup = get(SqChain_in, "Nup", 4)
    Ndn = get(SqChain_in, "Ndn", 4)
    t = get(SqChain_in, "t", 0.001)
    bias = get(SqChain_in, "bias", 0.0)
    
    
    run_SqChain(L, Nup, Ndn, t; U=U, bias=bias)

    # dyna_occ()
    # dyna_EE()
    # dyna_dptcurrent()

end 