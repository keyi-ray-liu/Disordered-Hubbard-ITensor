
function run_NF(L, Nup, Ndn, t;  U=4.0, bias=0.0, kwargs...)


    sys = set_NF_square(; L=L, Nup=Nup, Ndn=Ndn, t=t, U=U, bias=bias, kwargs...)

    simulation = set_Static(; sweepcnt=100, sweepdim = 300, kwargs...)
    ψ = gen_state(sys)

    run_static_simulation(sys, simulation, ψ)

end 



function NF_wrapper()

    NF_in = load_JSON( pwd() * "/NFpara.json")

    U = get(NF_in, "U", 4.0)
    L = get(NF_in, "L", 3)
    Nup = get(NF_in, "Nup", 4)
    Ndn = get(NF_in, "Ndn", 4)
    t = get(NF_in, "t", 0.001)
    bias = get(NF_in, "bias", 0.0)
    
    
    run_NF(L, Nup, Ndn, t; U=U, bias=bias)

    # dyna_occ()
    # dyna_EE()
    # dyna_dptcurrent()

end 