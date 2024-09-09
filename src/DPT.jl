"""get mixed basis reservoir parameters, the energies returned are at 0 bias, however the order is done at finite bias"""


gen_obs(mixed, includeU) = [dyna_EE, dyna_occ, (mixed && includeU) ? dyna_dptcurrent_mix : dyna_dptcurrent, dyna_corr, dyna_SRDM]


function continue_from_last()

    if !isempty(Glob.glob( "$LASTSTSTR.h5", getworkdir()))
        @show "non empty"
    else
        @show "empty"
    end 
end 


"""worker function that runs DPT calculations"""
function run_DPT(U, L, R, t_switch::Float64; bias_L = bias_LR/2, bias_R  = - bias_LR/2, τ=0.125, mixed=false, save_every=false,  ddposition="R", graph=false, avg=false, switchinterval ::Int=1,  kwargs...)


    eqinit_str = "EqInit"

    includeU = get(kwargs, :includeU, true)
    couple_range = get(kwargs, :couple_range, 2)
    ordering = get(kwargs, :ordering, "SORTED")
    t_doubledot = get(kwargs, :t_doubledot, 0.125 )

    # we first run a calculation with no bias on the LR, 

    energies, ks, LR = gen_mixed(mixed, L, R, bias_L, bias_R; ordering=ordering, includeU=includeU, couple_range = couple_range)
    obs = gen_obs(mixed, includeU)
    
    if L > 140
        filter!( e->e ∉ [dyna_corr], obs)
    end 

    if get(kwargs, :TEdim, 64) > 256 || L > 40
        save_every=false
    end 

    @show save_every
    @show obs

    eq = DPT_setter(mixed, avg; U=U, L=L, R=R, bias_doubledot=DPT_INIT_BIAS, t_doubledot=0.0, energies=energies, ks=ks, LR=LR, includeU=includeU, couple_range=couple_range, ddposition=ddposition, graph=graph)

    # if not present, we calculate the initial state

    if !check_ψ(eqinit_str)
        
        Static = set_Static(; output=eqinit_str, sweepdim=get(kwargs, :TEdim, 64) , kwargs...)

        # GS calculation
        ψ = gen_state(eq)
        ψ0 =  run_static_simulation(eq, Static, ψ; message = "Init")[1]
    else
        ψ0 = load_ψ(eqinit_str)
    end 

    # Stage1, no bias, GS, start negative time

    Stage1 = set_Dynamic(;τ=τ, start= - 1.0 + τ, fin=0, kwargs...)
    ψ1 = run_dynamic_simulation(eq, Stage1, ψ0; message="Stage1", save_every=save_every, obs=obs)

    # now we switch on the bias in L/R, 0 time

    noneq = DPT_setter(mixed, avg; U=U, L=L, R=R, t_doubledot=0.0, bias_L=bias_L, bias_R=bias_R, energies=energies, ks=ks, LR=LR, includeU=includeU, couple_range=couple_range, ddposition=ddposition, graph=graph)

    Stage2 = set_Dynamic(;τ=τ, start=τ, fin=t_switch , kwargs...)
    ψ2 = run_dynamic_simulation(noneq, Stage2, ψ1; message="Stage2", save_every=save_every, obs=obs, init_obs=false)

    # we then switch on the tunneling b/w drain_offset

    # if we have finite switch time in one timestep
    τ0 = τ/switchinterval
    for t in τ0:τ0:τ

        @show tempt_doubledot =  t_doubledot* t/τ 

        noneqtuninterval = DPT_setter(mixed, avg; U=U, L=L, R=R, t_doubledot=tempt_doubledot, bias_L=bias_L, bias_R=bias_R, energies=energies, ks=ks, LR=LR, includeU=includeU, couple_range=couple_range, ddposition=ddposition, graph=graph)

        Stage3interval = set_Dynamic(;τ=τ0, start=t_switch + t, fin=t_switch + t, kwargs...)

        ψ2 = run_dynamic_simulation(noneqtuninterval, Stage3interval, ψ2; message="Stage3interval",  save_every=save_every, obs=obs, init_obs=false)

    end 

    noneqtun = DPT_setter(mixed, avg; U=U, L=L, R=R, bias_L=bias_L, t_doubledot=t_doubledot, bias_R=bias_R, energies=energies, ks=ks, LR=LR, includeU=includeU, couple_range=couple_range, ddposition=ddposition, graph=graph)

    Stage3 = set_Dynamic(;τ=τ, start=t_switch + 2 * τ, fin=t_switch * 2, kwargs...)

    #ψ = load_ψ(t_switch)

    _ = run_dynamic_simulation(noneqtun, Stage3, ψ2; message="Stage3",  save_every=save_every, obs=obs, init_obs=false)

    return nothing
end 



function DPT_wrapper()

    

    dpt_in = load_JSON( pwd() * "/dptpara.json")

    U = get(dpt_in, "U", 2.0)
    L = get(dpt_in, "L", 16)
    R = get(dpt_in, "R", 16)
    t_switch = Float64(get(dpt_in, "tswitch", 5.0))
    τ = get(dpt_in, "timestep", 0.125)
    TEdim = get(dpt_in, "TEdim", 64)
    bias_LR = get(dpt_in, "bias_LR", BIAS_LR)
    mixed = get(dpt_in, "mixed", false)
    ordering = get(dpt_in, "ordering", "SORTED")
    includeU = get(dpt_in, "includeU", true)
    ddposition = get(dpt_in, "ddposition", "R")
    avg = get(dpt_in, "avg", false)
    switchinterval = get(dpt_in, "switchinterval", 20)
    sweepcnt = get(dpt_in, "sweepcnt", 60)
    t_doubledot = get(dpt_in, "tdoubledot", 0.125)

    if DISABLE_BLAS && TEdim > 256
        @show ITensors.enable_threaded_blocksparse()
        @show ITensors.enable_threaded_blocksparse()
    else
        @show ITensors.disable_threaded_blocksparse()
        @show ITensors.disable_threaded_blocksparse()
    end 

    run_DPT(U, L, R, t_switch; τ=τ, TEdim=TEdim, bias_L = bias_LR/2, bias_R  = -bias_LR/2, mixed=mixed, ordering=ordering, includeU=includeU, ddposition=ddposition, avg=avg, switchinterval=switchinterval, sweepcnt=sweepcnt, t_doubledot=t_doubledot)

    # dyna_occ()
    # dyna_EE()

    # if mixed
    #     dyna_dptcurrent_mix()
    # else
    #     dyna_dptcurrent()
    # end 

end 


function plot_mix()

    #sys = set_DPT_mixed(; L =8, R=8, graph=true )
    # sys = set_DPT(; L =6, R=6, graph=true )

    # ψ = gen_state(sys)
    L = 128
    energy, _, _ = gen_mixed(L, L, 0.25, -0.25; ordering="LRSORTED")
    energy, _, _ = gen_mixed(L, L, 0.25, -0.25; ordering="SORTED")
    energy, _, _ = gen_mixed(L, L, 0.25, -0.25; ordering="KMATCHED")
    energy, _, _ = gen_mixed(L, L, 0.25, -0.25; ordering="ABSSORTED")

    # open(getworkdir() * "testenergy", "w") do io
    #     writedlm(io, energy)
    # end 

    return nothing
end


function DPT_corr()

    dpt_in = load_JSON( pwd() * "/dptpara.json")
    avg = get(dpt_in, "avg", false)

    sys = DPT_setter(true, avg )
    dyna_corr(; sys=sys)
    dyna_SRDM()

end 