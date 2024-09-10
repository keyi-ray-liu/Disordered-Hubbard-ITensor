
"""worker function that runs DPT calculations"""
function run_LSR(L, R, fin; bias_L = BIAS_LR, bias_R = -BIAS_LR, τ=0.125, kwargs...)

    # we first run a calculation with no bias on the LR, 
    sys = set_LSR_SIAM(;L=L, R=R)
    eqinit_str = "EqInit"
    Static = set_Static(; output=eqinit_str)

    ψ = gen_state(sys)

    run_static_simulation(sys, Static, ψ)

    # now we switch on the bias in L/R
    noneq = set_LSR_SIAM(;L=L, R=R, bias_L=bias_L, bias_R=bias_R)

    Stage1 = set_Dynamic(;τ=τ, start=τ, fin=fin)
    ψ = load_ψ(eqinit_str)

    run_dynamic_simulation(noneq, Stage1, ψ)



end 


function LSR_SIAM_wrapper()

    lsr_in = load_JSON( pwd() * "/lsrpara.json")

    L = get(lsr_in, "L", 16)
    R = get(lsr_in, "R", 16)
    fin = get(lsr_in, "fin", 40.0)

    run_LSR(L, R, fin; τ=0.25)

    dyna_occ()
    #dyna_LSRcurrent()

end 



"""worker function that runs DPT calculations"""
function run_SD(L, R, fin; τ=0.125, bias_L=0.0, bias_R=0.0, kwargs...)
 
    obs= [dyna_EE, dyna_occ,  dyna_SRDM,
    #dyna_SDcurrent, dyna_corr,
    ]

    # we first run a calculation with no bias on the LR, 
    # sys = set_SD(;L=L, R=R, kwargs...)
    # eqinit_str = "EqInit"
    # Static = set_Static(; output=eqinit_str)

    #ψ = gen_state(sys)

    # run_static_simulation(sys, Static, ψ)

    # now we switch on the bias in L/R
    @show dyna = set_SD(;L=L, R=R, bias_L=bias_L, bias_R=bias_R, kwargs...)

    ψ = gen_state(dyna)
    Stage1 = set_Dynamic(;τ=τ, start=τ, fin=fin)
    #ψ = load_ψ(eqinit_str)

    _ = run_dynamic_simulation(dyna, Stage1, ψ; save_every=false, obs=obs)

    return nothing


end 


function SD_wrapper()

    sd_in = load_JSON( pwd() * "/sdpara.json")

    s_coupling = get(sd_in, "scoupling", -0.001)
    d_coupling = get(sd_in, "dcoupling", -0.001)
    λ_ne = get(sd_in, "int_ne", 2.0)
    λ_ee = get(sd_in, "int_ee", 2.0)
    L = get(sd_in, "L", 16)
    R = get(sd_in, "R", 16)
    Ns = get(sd_in, "Ns", 1)
    Nd = get(sd_in, "Nd", 0)
    Na = get(sd_in, "Na", 0)
    systype = get(sd_in, "systype", "Electron")
    fin = get(sd_in, "fin", 10.0)
    τ = get(sd_in, "timestep", 0.1)
    contact_scaling = get(sd_in, "contactscaling", 2.0)
    reservoir_type = get(sd_in, "reservoir_type", "spatial")

    run_SD(L, R, fin; τ=τ, Ns=Ns, Nd=Nd, Na=Na, s_coupling=s_coupling, d_coupling=d_coupling, 
    λ_ne = λ_ne, λ_ee = λ_ee, systype=systype, contact_scaling=contact_scaling, reservoir_type=reservoir_type)

end 