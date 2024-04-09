
"""worker function that runs DPT calculations"""
function run_DPT(U, L, R, t_switch::Float64; bias_L = BIAS_LR, bias_R = -BIAS_LR, τ=0.125, kwargs...)

    eqinit_str = "EqInit"
    # we first run a calculation with no bias on the LR, 
    eq = set_DPT(;U=U, L=L, R=R, bias_doubledot=DPT_INIT_BIAS, t_doubledot=0.0)

    # we make sure the sweepdim and TEdim match
    Static = set_Static(; output=eqinit_str, sweepdim=get(kwargs, :TEdim, 64), sweepcnt=100)

    ψ = gen_state(eq)

    run_static_simulation(eq, Static, ψ)

    # Stage1, no bias, GS, start negative time

    Stage1 = set_Dynamic(;τ=τ, start= - t_switch + τ, fin=t_switch, kwargs...)

    ψ = load_ψ(eqinit_str)

    run_dynamic_simulation(eq, Stage1, ψ)

    # now we switch on the bias in L/R, 0 time
    noneq = set_DPT(;U=U, L=L, R=R, t_doubledot=0.0,bias_L=bias_L, bias_R=bias_R)

    Stage2 = set_Dynamic(;τ=τ, start=τ, fin=t_switch , kwargs...)

    ψ = load_ψ(t_switch)

    run_dynamic_simulation(noneq, Stage2, ψ)

    # we then switch on the tunneling b/w drain_offset

    noneqtun = set_DPT(;U=U, L=L, R=R, bias_L=bias_L, bias_R=bias_R)

    Stage3 = set_Dynamic(;τ=τ, start=t_switch +τ, fin=t_switch * 2, kwargs...)

    ψ = load_ψ(t_switch)

    run_dynamic_simulation(noneqtun, Stage3, ψ)


end 


function DPT_wrapper()

    dpt_in = load_JSON( pwd() * "/dptpara.json")

    U = get(dpt_in, "U", 2.0)
    L = get(dpt_in, "L", 16)
    R = get(dpt_in, "R", 16)
    t_switch = Float64(get(dpt_in, "tswitch", 5.0))
    τ = get(dpt_in, "timestep", 0.125)
    TEdim = get(dpt_in, "TEdim", 64)
    

    run_DPT(U, L, R, t_switch; τ=τ, TEdim=TEdim)

    dyna_occ()
    dyna_EE()
    dyna_dptcurrent()

end 