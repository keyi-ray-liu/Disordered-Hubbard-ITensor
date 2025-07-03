
function run_NF(L, Nup, Ndn, t;  U=4.0, bias=0.0, kwargs...)


    sys = NF_square(; L=L, Nup=Nup, Ndn=Ndn, t=t, U=U, bias=bias, kwargs...)

    simulation = StaticSimulation(; sweepcnt=100, sweepdim = 300, kwargs...)
    ψ = gen_state(sys)

    workflag = ""

    run_static_simulation(sys, simulation, ψ, Identity(), workflag)

end 



function NF_wrapper()

    NF_in = load_JSON( pwd() * "/NFpara.json")

    U = get(NF_in, "U", 4.0)
    L = get(NF_in, "L", 3)
    Nup = get(NF_in, "Nup", 3)
    Ndn = get(NF_in, "Ndn", 4)
    t = get(NF_in, "t", -1.0)
    bias = get(NF_in, "bias", 0.0)
    
    
    run_NF(L, Nup, Ndn, t; U=U, bias=bias)

    # dyna_occ()
    # dyna_EE()
    # dyna_dptcurrent()

end 