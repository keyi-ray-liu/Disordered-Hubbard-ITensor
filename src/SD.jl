
"""worker function that runs DPT calculations"""
function run_LSR(L, R, fin; bias_L = biasLR, bias_R = -biasLR, τ=0.125, kwargs...)

    # we first run a calculation with no bias on the LR, 
    sys = set_LSR_SIAM(;L=L, R=R)
    Static = set_Static(; output=EQINIT_STR)

    ψ = gen_state(sys)

    run_static_simulation(sys, Static, ψ)

    # now we switch on the bias in L/R
    noneq = set_LSR_SIAM(;L=L, R=R, bias_L=bias_L, bias_R=bias_R)

    Stage1 = set_Dynamic(;τ=τ, start=τ, fin=fin)
    ψ = load_ψ(EQINIT_STR)

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


function run_SD(::ProductStateDriver, energies, ks, LR,  fin, τ, obs;  biasA=0.0, biasS=0.0, biasD=0.0, kwargs... )
    
    @show init = nothing
    @show sys = set_SD(; biasS = biasS, biasA = biasA, biasD=biasD, energies = energies, ks =ks, LR=LR, kwargs...)

    run_gs_dyna(τ, fin, init, sys, obs; kwargs...)
end 

function run_SD(::BiasReleaseDriver, energies, ks, LR, fin, τ, obs;  biasA=0.0, biasS=0.0, biasD=0.0, biasAinit=0.0, kwargs... )
    
    @show init = set_SD(; biasS = biasS, biasA = biasAinit, biasD = biasD, energies = energies, ks =ks, LR=LR, s_coupling = 0.0, d_coupling = 0.0, kwargs...)
    @show sys = set_SD(; biasS = 0.0, biasA = biasA, biasD=0.0, energies = energies, ks =ks, LR=LR, kwargs...)

    run_gs_dyna(τ, fin, init, sys, obs; kwargs...)
end 


function run_SD(::BiasGSDriver, energies, ks, LR, fin, τ, obs;  biasA=0.0, biasS=0.0, biasD=0.0, biasAinit=0.0, kwargs... )
    @assert biasS != biasD != 0

    @show init = set_SD(; biasS = 0, biasA = biasAinit, biasD = 0, energies = energies, ks =ks, LR=LR,  kwargs...)
    @show sys = set_SD(; biasS = biasS, biasA = biasA, biasD=biasD, energies = energies, ks =ks, LR=LR, kwargs...)

    run_gs_dyna(τ, fin, init, sys, obs; kwargs...)
end 


"""worker function that runs SD calculations"""
function run_SD(fin; τ=0.125, biasS=0.0, biasA=0.0, biasD=0.0, biasAinit = 500.0, mode="productstate",  kwargs...)
 
    obs= [dyna_EE, dyna_occ, dyna_SDcurrent,
    #dyna_SRDM, 
    ]

    # now we switch on the bias in L/R
    energies, ks, LR = gen_mixed( get(kwargs, :reservoir_type, "spatial")=="mixed", get(kwargs, :Ls, 4), get(kwargs, :Ld, 4), biasS, biasD; couple_range=0 )

    # else
    if mode == "productstate"
        modedriver = ProductStateDriver()

    elseif mode == "BiasRelease"
        modedriver = BiasReleaseDriver()

    elseif mode == "BiasGS"
        modedriver = BiasGSDriver()

    else
        error("Unrecognized Drive mode")
    end 

    run_SD(modedriver, energies, ks, LR, fin, τ, obs;  biasA=biasA, biasS=biasS, biasD=biasD, biasAinit=biasAinit, corr_cutoff=4.0, kwargs... )

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
    biasAinit = get(sd_in, "biasAinit", 500)
    systype = get(sd_in, "systype", "Electron")
    fin = get(sd_in, "fin", 10.0)
    τ = get(sd_in, "timestep", 0.1)
    contact_scaling = get(sd_in, "contactscaling", 2.0)
    reservoir_type = get(sd_in, "reservoirtype", "spatial")
    U = get(sd_in, "U", 4.0)
    config = get(sd_in, "config", "3x3")
    manualmixprod = get(sd_in, "manualmixprod", false)
    mode = get(sd_in, "mode", "productstate")

    run_SD(fin; τ=τ,  s_coupling=s_coupling, d_coupling=d_coupling, Ls=Ls, Ld=Ld, Ns=Ns, Na = Na, Nd=Nd, λ_ne = λ_ne, λ_ee = λ_ee, systype=systype, TEdim = TEdim, contact_scaling=contact_scaling, U=U, biasS = biasS, biasA = biasA, biasD = biasD, reservoir_type=reservoir_type, manualmixprod=manualmixprod, mode=mode, config=config, biasAinit = biasAinit)

end 