"""get mixed basis reservoir parameters, the energies returned are at 0 bias, however the order is done at finite bias"""
function gen_mixed(L, R, bias_L, bias_R; ordering="SORTED", includeU=true, couple_range=2)
    @info "Set mixed basis, $ordering"
    unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

    if !includeU
        L -= couple_range
        R -= couple_range
    end 

    # HACK 
    # we test individual version


    L_val = [ (2 *  cos( k * pi / (L + 1) )+ bias_L, k, 1) for k in 1:L] 
    R_val = [ (2 *  cos( k * pi / (R + 1) ) + bias_R, k, -1) for k in 1:R] 


    if ordering == "RNG"

        seed = abs(rand(Int))
        open(getworkdir * "randommixseed", "w") do io
            writedlm(io, seed)
        end 
        result = shuffle( StableRNG(seed), vcat(L_val, R_val))
        
    elseif ordering == "SORTED"
        result = sort( vcat(L_val, R_val), rev=true)

    elseif ordering == "MATCH"
        result = vcat(sort(L_val, by= e -> abs(e[1]), rev=true), sort(R_val, by= e -> abs(e[1])))

    else
        error("Unknown ordering")
    end

    energies, ks, LR= unzip(result)
    writedlm(getworkdir() * "mixenergies", energies)

    energies -= LR * bias_L

    writedlm(getworkdir() * "ks", ks)
    writedlm(getworkdir() * "LR", LR)

    return energies, ks, LR

end 



"""worker function that runs DPT calculations"""
function run_DPT(U, L, R, t_switch::Float64; bias_L = bias_LR/2, bias_R  = - bias_LR/2, τ=0.125, mixed=false, save_every=false,  ddposition="R", graph=false, kwargs...)


    eqinit_str = "EqInit"

    includeU = get(kwargs, :includeU, true)
    couple_range = get(kwargs, :couple_range, 2)
    ordering = get(kwargs, :ordering, "SORTED")

    # we first run a calculation with no bias on the LR, 

    if mixed
        energies, ks, LR = gen_mixed(L, R, bias_L, bias_R; ordering=ordering, includeU=includeU, couple_range = couple_range)
        eq = set_DPT_mixed(; U=U, L=L, R=R, bias_doubledot=DPT_INIT_BIAS, t_doubledot=0.0, energies=energies, ks=ks, LR=LR, includeU=includeU, couple_range=couple_range, ddposition=ddposition, graph=graph)

        current_obs = includeU ? dyna_dptcurrent_mix : dyna_dptcurrent
        obs = [dyna_EE, dyna_occ, current_obs, dyna_corr]
    else
        eq = set_DPT(;U=U, L=L, R=R, bias_doubledot=DPT_INIT_BIAS, t_doubledot=0.0, ddposition=ddposition, graph=graph)
        obs = [dyna_EE, dyna_occ, dyna_dptcurrent, dyna_corr]
    end 

    if get_systotal(eq) > 80
        filter!( e->e ∉ [dyna_corr], obs)
    end 

    @show obs

    # if not present, we calculate the initial state

    if !check_ψ(eqinit_str)
        
        Static = set_Static(; output=eqinit_str, sweepdim=get(kwargs, :TEdim, 64) , sweepcnt=80)

        # GS calculation
        ψ = gen_state(eq)
        ψ0 =  run_static_simulation(eq, Static, ψ; message = "Init")[1]
    else
        ψ0 = load_ψ(eqinit_str)
    end 

    # Stage1, no bias, GS, start negative time

    Stage1 = set_Dynamic(;τ=τ, start= - 4.0 + τ, fin=0, kwargs...)
    ψ1 = run_dynamic_simulation(eq, Stage1, ψ0; message="Stage1", save_every=save_every, obs=obs)

    # now we switch on the bias in L/R, 0 time

    if mixed
        noneq = set_DPT_mixed(;U=U, L=L, R=R, t_doubledot=0.0, bias_L=bias_L, bias_R=bias_R, energies=energies, ks=ks, LR=LR, includeU=includeU, couple_range=couple_range, ddposition=ddposition, graph=graph)
    else
        noneq = set_DPT(;U=U, L=L, R=R, t_doubledot=0.0,bias_L=bias_L, bias_R=bias_R, ddposition=ddposition, graph=graph)
    end 

    Stage2 = set_Dynamic(;τ=τ, start=τ, fin=t_switch , kwargs...)
    #ψ = load_ψ(0.0)

    ψ2 = run_dynamic_simulation(noneq, Stage2, ψ1; message="Stage2", save_every=save_every, obs=obs)

    # we then switch on the tunneling b/w drain_offset

    if mixed
        noneqtun = set_DPT_mixed(;U=U, L=L, R=R, bias_L=bias_L, bias_R=bias_R, energies=energies, ks=ks, LR=LR, includeU=includeU, couple_range=couple_range, ddposition=ddposition, graph=graph)
    else
        noneqtun = set_DPT(;U=U, L=L, R=R, bias_L=bias_L, bias_R=bias_R, ddposition=ddposition, graph=graph)
    end 

    Stage3 = set_Dynamic(;τ=τ, start=t_switch +τ, fin=t_switch * 2, kwargs...)

    #ψ = load_ψ(t_switch)

     _ = run_dynamic_simulation(noneqtun, Stage3, ψ2; message="Stage3",  save_every=save_every, obs=obs)

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

    run_DPT(U, L, R, t_switch; τ=τ, TEdim=TEdim, bias_L = bias_LR/2, bias_R  = -bias_LR/2, mixed=mixed, ordering=ordering, includeU=includeU, ddposition=ddposition)

    # dyna_occ()
    # dyna_EE()

    # if mixed
    #     dyna_dptcurrent_mix()
    # else
    #     dyna_dptcurrent()
    # end 

end 


function DPT_graph_test()

    #sys = set_DPT_mixed(; L =8, R=8, graph=true )
    # sys = set_DPT(; L =6, R=6, graph=true )

    # ψ = gen_state(sys)
    energy, _, _ = gen_mixed(128, 128, 0.25, -0.25; ordering="MATCH")

    open(getworkdir() * "testenergy", "w") do io
        writedlm(io, energy)
    end 


end 