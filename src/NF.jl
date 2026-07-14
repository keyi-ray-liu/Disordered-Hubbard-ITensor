
function run_NF(L, Nup, Ndn, t;  U=4.0, bias=0.0, flux=0.0, kwargs...)


    sys = NF_square(; L=L, Nup=Nup, Ndn=Ndn, t=t, U=U, bias=bias, flux=flux, kwargs...)

    simulation = StaticSimulation(ex=1; sweepcnt=30, sweepdim = 300, kwargs...)
    ψ = gen_state(sys)
    # ψ = load_ψ("wf.h5")


    run_static_simulation(sys, simulation, ψ, Identity())

end 



function NF_wrapper()

    NF_in = load_JSON( pwd() * "/NFpara.json")

    U = get(NF_in, "U", 4.0)
    L = get(NF_in, "L", 3)
    Nup = get(NF_in, "Nup", 4)
    Ndn = get(NF_in, "Ndn", 4)
    t = get(NF_in, "t", 0.001)
    bias = get(NF_in, "bias", 0.0)
    bool_rand = get(NF_in, "random_bool", false)
    flux = get(NF_in, "flux", 0.0)
    
    if bool_rand 
        random_onsite = (2 .* rand(L^2) .- 1) .* bias
        workdir = getworkdir()
        writedlm(workdir * "on_bias", random_onsite)
    end

    run_NF(L, Nup, Ndn, t; U=U, bias=bias, flux=flux)

    # dyna_occ()
    # dyna_EE()
    # dyna_dptcurrent()

end 