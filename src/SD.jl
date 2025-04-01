
"""worker function that runs DPT calculations"""
function run_LSR(L, R, fin; bias_L = biasLR, bias_R = -biasLR, τ=0.125, kwargs...)

    # we first run a calculation with no bias on the LR, 
    sys = set_LSR_SIAM(;L=L, R=R)
    Static = StaticSimulation(; output=EQINIT_STR)

    ψ = gen_state(sys)

    run_static_simulation(sys, Static, ψ, Identity())

    # now we switch on the bias in L/R
    noneq = set_LSR_SIAM(;L=L, R=R, bias_L=bias_L, bias_R=bias_R)

    Stage1 = DynamicSimulation(;τ=τ, start=τ, fin=fin)
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


function run_SD(::ProductStateDriver, timecontrol::TimeControl, energies, ks, LR, obs;  biasA=0.0, biasS=0.0, biasD=0.0, ω=-1.0, kwargs... )
    
    @show init = nothing
    @show sys = SD_array(; biasS = biasS, biasA = biasA, biasD=biasD, energies = energies, ks =ks, LR=LR, ω = ω,  kwargs...)

    run_gs_dyna(timecontrol, init, sys, obs; kwargs...)
end 

function run_SD(::ProdReservoirDriver, timecontrol::TimeControl, energies, ks, LR, obs;  Ns = 0, Nd = 0, reservoir_type = "mixed", ω = -1.0, kwargs... )
    
    if reservoir_type == "mixed"
        init = SD_array(;  energies = energies, ks =ks, LR=LR,  kwargs..., biasS = 1e16, biasA = 0.0, biasD = 1e16, Ns = 0, Nd = 0, ω = ω,  s_coupling = 0.0, d_coupling = 0.0, reservoir_type = "mixed")
        sys = SD_array(; energies = energies, ks =ks, LR=LR, kwargs..., Ns =Ns, Nd = Nd, biasS = 0.0, biasD=0.0, ω = ω,  reservoir_type = "mixed")

        process = LoadSource(Ns)
        run_gs_dyna(timecontrol, init, sys, obs; process = process, kwargs...)
    else
        init = SD_array(;  energies = energies, ks =ks, LR=LR,  kwargs..., biasS = -1e3, biasA = 0.0, biasD = 1e3, Ns = Ns, Nd = Nd, s_coupling = 0.0, d_coupling = 0.0, reservoir_type = "spatial", ω = ω )
        sys = SD_array(; energies = energies, ks =ks, LR=LR, kwargs..., biasS = 0.0, biasD=0.0, reservoir_type = "spatial", ω = ω)

        run_gs_dyna(timecontrol, init, sys, obs;  kwargs...)

    end 
end 



function run_SD(::BiasReverseGS, timecontrol::TimeControl, energies, ks, LR, obs;  kwargs... )
    
    init = SD_array(;  energies = energies, ks =ks, LR=LR,  kwargs..., )
    sys = SD_array(; energies = energies, ks =ks, LR=LR, kwargs...,  biasS = 0.0, biasD=0.0)

    run_gs_dyna(timecontrol, init, sys, obs;  kwargs...)

end 

function run_SD(::BiasGSDriver, timecontrol::TimeControl, energies, ks, LR, obs;  biasA=0.0, biasS=0.0, biasD=0.0, initbiasA=0.0, ω = -1.0, kwargs... )
    @assert biasS != biasD != 0

    @show init = SD_array(; biasS = 0, biasA = initbiasA, biasD = 0, energies = energies, ks =ks, LR=LR, ω = ω,  kwargs...)
    @show sys = SD_array(; biasS = biasS, biasA = biasA, biasD=biasD, energies = energies, ks =ks, LR=LR, ω = ω,  kwargs...)

    run_gs_dyna(timecontrol, init, sys, obs; kwargs...)
end 


"""worker function that runs SD calculations"""
function run_SD(; biasS=0.0, biasA=0.0, biasD=0.0, initbiasA = 500.0, mode="productstate", ω = -1.0, kwargs...)
 
    obs= [dyna_EE, dyna_occ, dyna_SDcurrent,
    #dyna_SRDM, 
    #dyna_product_state_overlap
    ]

    # now we switch on the bias in L/R
    
    # else
    if mode == "productstate"
        energies, ks, LR = gen_mixed( get(kwargs, :reservoir_type, "spatial")=="mixed"; L = get(kwargs, :Ls, 4), R = get(kwargs, :Ld, 4), bias_L = biasS, bias_R = biasD, couple_range=0, ω = ω )
        modedriver = ProductStateDriver()

    elseif mode == "ProdReservoir"
        energies, ks, LR = gen_mixed( get(kwargs, :reservoir_type, "spatial")=="mixed"; L = get(kwargs, :Ls, 4), R = get(kwargs, :Ld, 4), bias_L = 0.0, bias_R = 0.0, couple_range=0, ω = ω )
        modedriver = ProdReservoirDriver()

    elseif mode == "BiasReverseGS"
        energies, ks, LR = gen_mixed( get(kwargs, :reservoir_type, "spatial")=="mixed"; L = get(kwargs, :Ls, 4), R = get(kwargs, :Ld, 4), bias_L = 0.0, bias_R = 0.0, couple_range=0, ω = ω )
        modedriver = BiasReverseGS()

    elseif mode == "BiasGS"

        energies, ks, LR = gen_mixed( get(kwargs, :reservoir_type, "spatial")=="mixed"; L = get(kwargs, :Ls, 4), R = get(kwargs, :Ld, 4), bias_L = biasS, bias_R = biasD, couple_range=0, ω = ω)
        modedriver = BiasGSDriver()

    else
        error("Unrecognized Drive mode")
    end 

    
    timecontrol = get_time_control()

    run_SD(modedriver, timecontrol, energies, ks, LR, obs;  biasA=biasA, biasS=biasS, biasD=biasD, initbiasA=initbiasA, corr_cutoff=0.0, ω=ω, kwargs... )

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
    initbiasA = get(sd_in, "initbiasA", 500)
    systype = get(sd_in, "systype", "Electron")
    contact_scaling = get(sd_in, "contactscaling", 2.0)
    reservoir_type = get(sd_in, "reservoirtype", "spatial")
    U = get(sd_in, "U", 4.0)
    config = get(sd_in, "config", "3x3")
    manualmixprod = get(sd_in, "manualmixprod", false)
    mode = get(sd_in, "mode", "productstate")
    ω = get(sd_in, "reservoirspacing", -1.0)
    sweepcnt = get(sd_in, "sweepcnt", 200)

    run_SD(; s_coupling=s_coupling, d_coupling=d_coupling, Ls=Ls, Ld=Ld, Ns=Ns, Na = Na, Nd=Nd, λ_ne = λ_ne, λ_ee = λ_ee, systype=systype, TEdim = TEdim, contact_scaling=contact_scaling, U=U, biasS = biasS, biasA = biasA, biasD = biasD, reservoir_type=reservoir_type, manualmixprod=manualmixprod, mode=mode, config=config, initbiasA = initbiasA, ω=ω, sweepcnt = sweepcnt)

end 