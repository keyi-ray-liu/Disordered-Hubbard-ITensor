function run_rectangle(Lx, Ly, Nup, Ndn, t;  U=4.0, V=0.5, bias=0.0, flux=flux, kwargs...)


    sys = NF_rect(; Lx=Lx, Ly=Ly, Nup=Nup, Ndn=Ndn, t=t, U=U, V=V, bias=bias, flux=flux, kwargs...)

    simulation = StaticSimulation(ex=1; sweepcnt=100, sweepdim = 100, kwargs...)
    # ψ = gen_state(sys)
    ψ = load_ψ("wf.h5")

    run_static_simulation(sys, simulation, ψ, Identity())

end 


function rect_wrapper()

    rect_in = load_JSON( pwd() * "/rect_para.json")

    U = get(rect_in, "U", 4.0)
    V = get(rect_in, "V", 0.5)
    flux = get(rect_in, "flux", 0.0)
    Lx = get(rect_in, "Lx", 3)
    Ly = get(rect_in, "Ly", 4)
    Nup = get(rect_in, "Nup", 4)
    Ndn = get(rect_in, "Ndn", 4)
    t = get(rect_in, "t", 0.001)
    bias = get(rect_in, "bias", 0.0) 
    bool_rand = get(rect_in, "random_bool", false) 
    
    if bool_rand 
        random_onsite = (2 .* rand(Lx*Ly) .- 1) .* bias
        workdir = getworkdir()
        writedlm(workdir * "on_bias", random_onsite)
    end

    run_rectangle(Lx, Ly, Nup, Ndn, t; U=U, V=V, bias=bias, flux=flux)

    # dyna_occ()
    # dyna_EE()
    # dyna_dptcurrent()

end 