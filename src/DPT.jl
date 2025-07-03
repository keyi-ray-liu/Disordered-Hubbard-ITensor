"""get mixed basis Reservoir parameters, the energies returned are at 0 bias, however the order is done at finite bias"""


gen_obs(mixed, QPCmixed) = [dyna_EE, dyna_occ, dyna_coherence, (mixed && QPCmixed) ? dyna_dptcurrent_mix : dyna_dptcurrent, dyna_corr, dyna_SRDM]








"""worker function that runs DPT calculations"""
function run_DPT_many_body(U, L, R,  t_fin :: Float64; tswitch = 0.0, bias_L = BIASLR/2, bias_R  = - BIASLR/2, τ=0.25, mixed=false,  ddposition="R", graph=false, avg=false,  n1init = 0.0, QPC = 0.0, workflag = "", vs = 0.25, mode = "disconnectDD", ordering = "SORTED", process = Identity(), sites = nothing , initdd = "LOWER", kwargs...)

    
    QPCmixed = get(kwargs, :QPCmixed, false)
    couple_range = get(kwargs, :couple_range, 2)


    # we first run a calculation with no bias on the LR, 

    energies, ks, LR = gen_mixed( get(kwargs, :reservoir_type, "mixed")=="mixed"; L = L, R = R,  ordering = ordering, bias_L = bias_L, bias_R = bias_R, couple_range=couple_range, ω = 1.0, workflag = workflag, QPCmixed = QPCmixed)
    obs = gen_obs(mixed, QPCmixed)

    if mode == "disconnectDD"

        μ1 = 0.0

        @info "Initial state with no correlation b/w dd: U = $(U)"

        Uinit = U
        vsinit = 1e-14

        DPT_INIT_BIAS = [-100, 0]

    elseif mode == "connectDD"

        μ1 = U * (n1init - 1/2)
        μC = U * (QPC - 2)

        @info "Initial state with correlation b/w dd: U = 0, μC = $(μC), μ1 = $(μ1)"
        Uinit = 0.0
        vsinit = vs

        DPT_INIT_BIAS = [μC, 0]

    elseif mode == "manuallower"

        μ1 = U * (n1init - 1/2)
        Uinit = 0.0

        vsinit = 1e-14
        DPT_INIT_BIAS = [-100, 0]

    elseif mode == "manualupper"


        μ1 = U * (n1init - 1/2)
        Uinit = 0.0

        vsinit = 1e-14
        DPT_INIT_BIAS = [0, -100]

    elseif mode == "manualempty"


        μ1 = U * (n1init - 1/2)
        Uinit =  0.0
        vsinit = 0.0
        DPT_INIT_BIAS = [+1e3, +1e3]

    elseif mode == "override"

        μ1 = 0.0
        Uinit = 0.0
        vsinit = 0.0
        DPT_INIT_BIAS = [0, -100]

        @assert typeof(process) == Supplywf

    else
        error("Unknown mode for DPT driving")
    end 



    # init we solve the decoupled system
    init = DPT_setter(mixed, avg; U=Uinit, L=L, R=R, bias_L = 0.0, bias_R = 0.0, bias_doubledot=DPT_INIT_BIAS, vs=vsinit, energies=energies, ks=ks, LR=LR, QPCmixed=QPCmixed, couple_range=couple_range, ddposition=ddposition, graph=graph, μ1 = μ1)

    if tswitch > 0

        stage2 = DPT_setter(mixed, avg; U=U, L=L, R=R, bias_L=bias_L, vs=0.0, bias_R=bias_R, bias_doubledot = [0.0, 0.0], energies=energies, ks=ks, LR=LR, QPCmixed=QPCmixed, couple_range=couple_range, ddposition=ddposition, graph=graph)
        stage3 = DPT_setter(mixed, avg; U=U, L=L, R=R, bias_L=bias_L, vs=vs, bias_R=bias_R, energies=energies, ks=ks, LR=LR, QPCmixed=QPCmixed, couple_range=couple_range, ddposition=ddposition, graph=graph)
        timecontrol = TwoStage(τ, τ, tswitch, t_fin)
        ψ = run_gs_dyna(timecontrol, init, stage2, stage3, obs; process = process, workflag = workflag, sites = sites, initdd = initdd, kwargs...)
    else
        sys = DPT_setter(mixed, avg; U=U, L=L, R=R, bias_L=bias_L, vs=vs, bias_R=bias_R, bias_doubledot = [0.0, 0.0], energies=energies, ks=ks, LR=LR, QPCmixed=QPCmixed, couple_range=couple_range, ddposition=ddposition, graph=graph)
        timecontrol = OneStage( τ, t_fin)
        ψ =  run_gs_dyna(timecontrol, init, sys, obs; process = process, workflag = workflag, sites = sites, initdd = initdd, kwargs...)
    end 
    
    

    return ψ
end 




function DPT_wrapper()

    

    dpt_in = load_JSON( pwd() * "/dptpara.json")

    U = get(dpt_in, "U", 2.0)
    L = get(dpt_in, "L", 16)
    R = get(dpt_in, "R", 16)
    tswitch = Float64(get(dpt_in, "tswitch", 5.0))
    t_fin = Float64(get(dpt_in, "tfin", 2 * tswitch))
    τ = get(dpt_in, "timestep", 0.125)
    TEdim = get(dpt_in, "TEdim", 64)
    biasLR = get(dpt_in, "biasLR", BIASLR)
    mixed = get(dpt_in, "mixed", false)
    ordering = get(dpt_in, "ordering", "SORTED")
    QPCmixed = get(dpt_in, "QPCmixed", true)
    ddposition = get(dpt_in, "ddposition", "R")
    avg = get(dpt_in, "avg", false)
    switchinterval = get(dpt_in, "switchinterval", 20)
    sweepcnt = get(dpt_in, "sweepcnt", 60)
    vs = get(dpt_in, "vs", 0.125)
    initdd = get(dpt_in, "initdd", "LOWER")
    mode = get(dpt_in, "mode", "disconnectDD")

    if DISABLE_BLAS && TEdim > 128  
        @show ITensors.enable_threaded_blocksparse()
        @show ITensors.enable_threaded_blocksparse()
    else
        @show ITensors.disable_threaded_blocksparse()
        @show ITensors.disable_threaded_blocksparse()
    end 

    #run_DPT(U, L, R, tswitch, t_fin; τ=τ, TEdim=TEdim, bias_L = biasLR/2, bias_R  = -biasLR/2, mixed=mixed, ordering=ordering, QPCmixed=QPCmixed, ddposition=ddposition, avg=avg, switchinterval=switchinterval, sweepcnt=sweepcnt, initdd = initdd, vs=vs)

    d1, QPC = get_mf(mode ; U = U, mu = biasLR, L = L, vs = vs )

    run_DPT_many_body(U, L, R,  t_fin; tswitch = tswitch, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  ddposition="R",  avg=avg,  n1init = d1, QPC = QPC, TEdim = TEdim, sweepcnt = sweepcnt, mode = mode, vs = vs, ordering = ordering)
    # dyna_occ()
    # dyna_EE()

    # if mixed
    #     dyna_dptcurrent_mix()
    # else
    #     dyna_dptcurrent()
    # end 

end 


function plot_mix()

    #sys = DPT_mixed(; L =8, R=8, graph=true )
    # sys = DPT(; L =6, R=6, graph=true )

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


function get_d1QPC(workflag, L)
    wf = open(getworkdir(workflag) * "occ", "r") 
    occ = readdlm(wf)
    close(wf)

    nd1 = sum(mean(occ[(end - 8):end , end - 1], dims = 1))
    nQPC = sum(mean(occ[(end - 8):end, (L - 1):(L + 2)], dims = 1))

    @show nd1, nQPC
    return nd1, nQPC
end 

function DPT_init_scan()



    dpt_in = load_JSON( pwd() * "/dptpara.json")
    U = get(dpt_in, "U", 0.1)
    L = get(dpt_in, "L", 34)
    R = get(dpt_in, "R", L)
    tswitch = Float64(get(dpt_in, "tswitch", 0.0))
    t_fin = Float64(get(dpt_in, "tfin", 34.0))
    τ = get(dpt_in, "timestep", 0.25)
    TEdim = get(dpt_in, "TEdim", 64)
    biasLR = get(dpt_in, "biasLR", 0.0)
    mixed = get(dpt_in, "mixed", true)
    ordering = get(dpt_in, "ordering", "SORTED")
    avg = get(dpt_in, "avg", false)
    sweepcnt = get(dpt_in, "sweepcnt", 30)
    vs = get(dpt_in, "vs", 0.25)


    αs = [0.5]
    #αs = 0.55:0.05:1.0

    for α in αs
        ψ = run_DPT_many_body(U, L, R,  0.0; tswitch = 0.0, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  ddposition="R",  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "manualempty", n1init = α, vs = vs, ordering = ordering, workflag = "zero_U$(U)_init$(α)", initdd = "EMPTY")
        
        ψ = two_site_rotate(ψ, α, 0)
        expval = expect(ψ, "N")
        corr = correlation_matrix(ψ, "Cdag", "C")

        @show expval[end - 1: end]
        @show corr[end - 1: end, end - 1: end]

        process = Supplywf(ψ)

        workflag = "_U$(U)_init$(α)"

        ψ2 = run_DPT_many_body(U, L, R,  t_fin; tswitch = tswitch, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  workflag = workflag, ddposition="R",  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "override", vs = vs, process = process)

        last = L + R + 2
        rdm = RDM2(ψ2, last - 1, last)

        # c1 = op("Cdag", siteinds(ψ2)[end - 1])
        # c2 = op("C", siteinds(ψ2)[end])

        workdir = getworkdir(workflag)
        h5open(workdir *"2RDM.h5", "w" ) do io
            write(io, "rho", rdm)
        end 

        
    end 


end 




function DPT_trend()


    dpt_in = load_JSON( pwd() * "/dptpara.json")
    U = get(dpt_in, "U", 0.1)
    L = get(dpt_in, "L", 34)
    R = get(dpt_in, "R", L)
    tswitch = Float64(get(dpt_in, "tswitch", 0.0))
    t_fin = Float64(get(dpt_in, "tfin", 34.0))
    τ = get(dpt_in, "timestep", 0.25)
    TEdim = get(dpt_in, "TEdim", 64)
    biasLR = get(dpt_in, "biasLR", 0.0)
    mixed = get(dpt_in, "mixed", true)
    ordering = get(dpt_in, "ordering", "SORTED")
    avg = get(dpt_in, "avg", false)
    sweepcnt = get(dpt_in, "sweepcnt", 30)
    vs = get(dpt_in, "vs", 0.25)


    d1s = []
    stepsize = 0.2
    Us = [ U - stepsize, U]
    # establish initial guess

    for U in Us

        workflag = "_init_U" * string(trunc(U, sigdigits=5)) 
        run_DPT_many_body(U, L, R,  t_fin; tswitch = tswitch, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  workflag = workflag, ddposition="R",  avg=avg, TEdim = TEdim, sweepcnt = sweepcnt, mode = "disconnectDD", vs = vs)

        nd1, _ = get_d1QPC(workflag, L)
        append!(d1s, sqrt(nd1))

    end 

    # @. linear(x, p) = p[1] * x + p[2]
    # p0 = [0.1, 0.5]

    # d1fit = curve_fit(linear, Us, d1s, p0)
    # @show coef(d1fit)

    # d1func(x) = linear(x, coef(d1fit))


    U = U - stepsize * 2
    while U > 0.01

        U = trunc(U, sigdigits=5)
        α = (max(1/sqrt(2), d1s[1] - (d1s[2] - d1s[1])))^2

        @show d1s, α
        workflag = "_U$(U)" 

        ψ = run_DPT_many_body(U, L, R,  0.0; tswitch = 0.0, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  ddposition="R",  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "manualempty", n1init = α, vs = vs, ordering = ordering, workflag = "zero_U$(U)_init$(α)", initdd = "EMPTY")
        
        ψ = two_site_rotate(ψ, α, 0)
        expval = expect(ψ, "N")
        corr = correlation_matrix(ψ, "Cdag", "C")

        @show expval[end - 1: end]
        @show corr[end - 1: end, end - 1: end]

        process = Supplywf(ψ)

        workflag = "_U$(U)_init$(α)"

        run_DPT_many_body(U, L, R,  t_fin; tswitch = tswitch, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  workflag = workflag, ddposition="R",  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "override", vs = vs, process = process)

        nd1, _ = get_d1QPC(workflag, L)
        d1s = [sqrt(nd1), d1s[1]]

        U -= stepsize
    end 

end 
