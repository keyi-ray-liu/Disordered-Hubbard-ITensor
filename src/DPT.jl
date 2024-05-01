"""get mixed basis reservoir parameters, the energies returned are at 0 bias, however the order is done at finite bias"""
function gen_mixed(L, R, bias_L, bias_R; random=false)

    unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))
    L_val = [ (2 *  cos( k * pi / (L + 1) )+ bias_L, k, 1) for k in 1:L] 
    R_val = [ (2 *  cos( k * pi / (R + 1) ) + bias_R, k, -1) for k in 1:R] 

    if random
        result = shuffle( vcat(L_val, R_val))
    else
        result = sort( vcat(L_val, R_val), rev=true)
    end

    energies, ks, LR= unzip(result)
    
    energies -= LR * bias_L

    writedlm(getworkdir() * "ks", ks)
    writedlm(getworkdir() * "LR", LR)

    return energies, ks, LR

end 



"""worker function that runs DPT calculations"""
function run_DPT(U, L, R, t_switch::Float64; bias_L = bias_LR/2, bias_R  = - bias_LR/2, τ=0.125, mixed=false, random=false, kwargs...)


    eqinit_str = "EqInit"
    # we first run a calculation with no bias on the LR, 

    if mixed
        energies, ks, LR = gen_mixed(L, R, bias_L, bias_R; random=random)
        eq = set_DPT_mixed(; U=U, L=L, R=R, bias_doubledot=DPT_INIT_BIAS, t_doubledot=0.0, energies=energies, ks=ks, LR=LR)
    else
        eq = set_DPT(;U=U, L=L, R=R, bias_doubledot=DPT_INIT_BIAS, t_doubledot=0.0)
    end 

    # we make sure the sweepdim and TEdim match

    if !check_ψ(eqinit_str)
        
        Static = set_Static(; output=eqinit_str, sweepdim=get(kwargs, :TEdim, 64) * 2, sweepcnt=100, cutoff=1E-11)

        # GS calculation
        ψ = gen_state(eq)
        run_static_simulation(eq, Static, ψ; message = "No init file found. Start init calculation")
    end 

    # Stage1, no bias, GS, start negative time

    Stage1 = set_Dynamic(;τ=τ, start= - 4.0 + τ, fin=0, kwargs...)

    ψ = load_ψ(eqinit_str)

    run_dynamic_simulation(eq, Stage1, ψ; message="Stage 1 begin")

    # now we switch on the bias in L/R, 0 time

    if mixed
        noneq = set_DPT_mixed(;U=U, L=L, R=R, t_doubledot=0.0, bias_L=bias_L, bias_R=bias_R, energies=energies, ks=ks, LR=LR)
    else
        noneq = set_DPT(;U=U, L=L, R=R, t_doubledot=0.0,bias_L=bias_L, bias_R=bias_R)
    end 

    Stage2 = set_Dynamic(;τ=τ, start=τ, fin=t_switch , kwargs...)

    ψ = load_ψ(0.0)

    run_dynamic_simulation(noneq, Stage2, ψ; message="Stage 2 begin")

    # we then switch on the tunneling b/w drain_offset

    if mixed
        noneqtun = set_DPT_mixed(;U=U, L=L, R=R, bias_L=bias_L, bias_R=bias_R, energies=energies, ks=ks, LR=LR)
    else
        noneqtun = set_DPT(;U=U, L=L, R=R, bias_L=bias_L, bias_R=bias_R)
    end 

    Stage3 = set_Dynamic(;τ=τ, start=t_switch +τ, fin=t_switch * 2, kwargs...)

    ψ = load_ψ(t_switch)

    run_dynamic_simulation(noneqtun, Stage3, ψ; message="Stage 3 begin")


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
    random = get(dpt_in, "random", false)
    

    run_DPT(U, L, R, t_switch; τ=τ, TEdim=TEdim, bias_L = bias_LR/2, bias_R  = -bias_LR/2, mixed=mixed, random=random)

    dyna_occ()
    dyna_EE()

    if mixed
        dyna_dptcurrent_mix()
    else
        dyna_dptcurrent()
    end 

end 