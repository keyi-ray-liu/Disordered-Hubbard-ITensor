
"""worker function that runs DPT calculations"""
function run_LSR(L, R, fin; bias_L = biasLR, bias_R = -biasLR, τ=0.125, kwargs...)

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



"""worker function that runs SD calculations"""
function run_SD(fin; τ=0.125, biasS=0.0, biasA=0.0, biasD=0.0,manualmixprod=false, mode="productstate",  kwargs...)
 
    obs= [dyna_EE, dyna_occ,  dyna_SRDM, dyna_corr
    #dyna_SDcurrent, dyna_corr,
    ]

    # we first run a calculation with no bias on the LR, 
    # sys = set_SD(;L=L, R=R, kwargs...)
    # eqinit_str = "EqInit"
    # Static = set_Static(; output=eqinit_str)

    #ψ = gen_state(sys)

    # run_static_simulation(sys, Static, ψ)

    # now we switch on the bias in L/R
    energies, ks, LR = gen_mixed( get(kwargs, :reservoir_type, "spatial")=="mixed", get(kwargs, :Ls, 4), get(kwargs, :Ld, 4), biasS, biasD; couple_range=0 )

    @show sys = set_SD(; biasA = biasA, biasS = biasS, biasD=biasD, energies = energies, ks =ks, LR =LR, kwargs...)

    ψ = gen_state(sys, manualmixprod=manualmixprod)

    
    # if ED
    #     run_exact_diagonalization(sys, ψ)

    # else
    if mode == "productstate"
        Stage1 = set_Dynamic(;τ=τ, start=τ, fin=fin, kwargs...)
        #ψ = load_ψ(eqinit_str)

        _ = run_dynamic_simulation(sys, Stage1, ψ; save_every=false, obs=obs)
    elseif mode == "leftGS"

        eqinit_str = "Eqinit"
        # we have a extremely strong bias on Left so we can load left
        init = set_SD(; biasS = -1000, biasA = 0, biasD = 0, energies = energies, ks =ks, LR=LR, kwargs...)

        Static = set_Static(; output=eqinit_str, sweepdim=get(kwargs, :TEdim, 64), ex=1, kwargs...)

        # GS calculation
        ψ0 =  run_static_simulation(init, Static, ψ; message = "Init")[1]

        Dynamic = set_Dynamic(;τ=τ, start=τ, fin=fin, kwargs...)
        _ = run_dynamic_simulation(sys, Dynamic, ψ0; save_every=false, obs=obs)

    elseif mode == "LRbias"

        @assert biasS != biasD != 0

        eqinit_str = "Eqinit"
        # we solve for the GS of the whole system at zero bias, we also bias the array so that nothing is occupied there
        init = set_SD(; biasS = 0, biasA = 500, biasD = 0, energies = energies, ks =ks, LR=LR, kwargs...)

        Static = set_Static(; output=eqinit_str, sweepdim=get(kwargs, :TEdim, 64), sweepcnt=80, ex=1, kwargs...)

        # GS calculation
        ψ0 =  run_static_simulation(init, Static, ψ; message = "Init")[1]

        Dynamic = set_Dynamic(;τ=τ, start=τ, fin=fin, kwargs...)
        _ = run_dynamic_simulation(sys, Dynamic, ψ0; save_every=false, obs=obs)

    end 


    # end 

    return nothing


end 


function SD_wrapper()

    sd_in = load_JSON( pwd() * "/sdpara.json")

    s_coupling = get(sd_in, "scoupling", -0.001)
    d_coupling = get(sd_in, "dcoupling", -0.001)
    λ_ne = get(sd_in, "intne", 0.0)
    λ_ee = get(sd_in, "intee", 0.0)
    TEdim = get(sd_in, "TEdim", 64)
    Ls = get(sd_in, "Ls", 1)
    Ld = get(sd_in, "Ld", Ls)
    Ns = get(sd_in, "Ns", 1)
    Nd = get(sd_in, "Nd", 0)
    Na = get(sd_in, "Na", 0)
    biasS = get(sd_in, "biasS", 0.0)
    biasA = get(sd_in, "biasA", 0.0)
    biasD = get(sd_in, "biasD", 0.0)
    systype = get(sd_in, "systype", "Electron")
    fin = get(sd_in, "fin", 10.0)
    τ = get(sd_in, "timestep", 0.1)
    contact_scaling = get(sd_in, "contactscaling", 2.0)
    reservoir_type = get(sd_in, "reservoir_type", "spatial")
    U = get(sd_in, "U", 4.0)
    manualmixprod = get(sd_in, "manualmixprod", false)
    mode = get(sd_in, "mode", "productstate")

    run_SD(fin; τ=τ,  s_coupling=s_coupling, d_coupling=d_coupling, Ls=Ls, Ld=Ld, Ns=Ns, Na = Na, Nd=Nd, 
    λ_ne = λ_ne, λ_ee = λ_ee, systype=systype, TEdim = TEdim, contact_scaling=contact_scaling, U=U, biasS = biasS, biasA = biasA, biasD = biasD, reservoir_type=reservoir_type, manualmixprod=manualmixprod, mode=mode)

end 