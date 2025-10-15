"""get mixed basis Reservoir parameters, the energies returned are at 0 bias, however the order is done at finite bias"""


gen_obs(mixed, QPCmixed) = [dyna_EE, dyna_occ, (mixed && QPCmixed) ? dyna_dptcurrent_mix : dyna_dptcurrent,
#, dyna_corr, dyna_SRDM 
# dyna_coherence,
dyna_SVD
]








"""worker function that runs DPT calculations"""
function run_DPT_many_body(U, L, R,  t_fin :: Float64; tswitch = 0.0, bias_L = BIASLR/2, bias_R  = - BIASLR/2, τ=0.25, mixed=false,  ddposition="R", graph=false, avg=false,  TLS = false, n1init = 0.0, QPC = 0.0, workflag = "", vs = 0.25, mode = "disconnectDD", ordering = "SORTED", process = Identity(), sites = nothing , initdd = "LOWER", n1penalty = nothing, kwargs...)

    
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

    @show μ1


    # init we solve the decoupled system
    init = DPT_setter(mixed, avg, TLS; U=Uinit, L=L, R=R, bias_L = 0.0, bias_R = 0.0, bias_doubledot=DPT_INIT_BIAS, vs=vsinit, energies=energies, ks=ks, LR=LR, QPCmixed=QPCmixed, couple_range=couple_range, ddposition=ddposition, graph=graph, μ1 = μ1)

    if tswitch > 0

        stage2 = DPT_setter(mixed, avg, TLS; U=U, L=L, R=R, bias_L=bias_L, vs=0.0, bias_R=bias_R, bias_doubledot = [0.0, 0.0], energies=energies, ks=ks, LR=LR, QPCmixed=QPCmixed, couple_range=couple_range, ddposition=ddposition, graph=graph, n1penalty = n1penalty)

        stage3 = DPT_setter(mixed, avg, TLS; U=U, L=L, R=R, bias_L=bias_L, vs=vs, bias_R=bias_R, bias_doubledot = [0.0, 0.0], energies=energies, ks=ks, LR=LR, QPCmixed=QPCmixed, couple_range=couple_range, ddposition=ddposition, graph=graph)

        timecontrol = TwoStage(τ, τ, tswitch, t_fin)

        ψ = run_gs_dyna(timecontrol, init, stage2, stage3, obs; process = process, workflag = workflag, sites = sites, initdd = initdd, kwargs...)
    else
        sys = DPT_setter(mixed, avg, TLS; U=U, L=L, R=R, bias_L=bias_L, vs=vs, bias_R=bias_R, bias_doubledot = [0.0, 0.0], energies=energies, ks=ks, LR=LR, QPCmixed=QPCmixed, couple_range=couple_range, ddposition=ddposition, graph=graph)

        timecontrol = OneStage( τ, t_fin)

        ψ =  run_gs_dyna(timecontrol, init, sys, obs; process = process, workflag = workflag, sites = sites, initdd = initdd, kwargs...)
    end 
    
    

    return ψ
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


function get_d1QPC(workflag, L, ddsite)
    wf = open(getworkdir(workflag) * "occ", "r") 
    occ = readdlm(wf)
    close(wf)

    nd1 = sum(mean(occ[(end - 16):end , ddsite], dims = 1))
    nQPC = sum(mean(occ[(end - 16):end, (L - 1):(L + 2)], dims = 1))

    @show workflag, nd1, nQPC
    return nd1, nQPC
end 

function DPT_wrapper(; dpt_in = nothing)


    if isnothing(dpt_in)
        dpt_in = load_JSON( pwd() * "/dptpara.json")
    end 
    
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
    repeat = get(dpt_in, "repeat", 1)
    TLS = get(dpt_in, "TLS", false)
    stagetype = get(dpt_in, "stagetype", "uniform")
    ddposition = get(dpt_in, "ddposition", "M")
    QN = get(dpt_in, "QN", true)
    #initlinkdim = get(dpt_in, "initlinkdim", 1)

    # the argument has higher priority
    α = get(dpt_in, "n1init", 1.0)

    if avg
        ddposition = "avg"
    end 

    for cnt in 1:repeat
        @info "reapet counter:", cnt
        @show α
        #workflag = "zero_repeat$(cnt)"
        workflag = "zero_repeat$(cnt)"

        ψ = run_DPT_many_body(U, L, R,  0.0; tswitch = 0.0, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  ddposition=ddposition,  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "manualempty", n1init = α, vs = vs, ordering = ordering, workflag = workflag, initdd = "EMPTY", TLS = TLS, QN = QN)
        
        #exppre = expect(ψ, "N")
        if ddposition == "R"
            ddsite = L + R + 1
        elseif ddposition == "M"
            ddsite = L + 1
        elseif ddposition == "avg"
            ddsite = L
        else
            ddsite = 1
        end 

        ψ = rotate_gate(ψ, α, 0, ddsite, TLS, avg)

        exppost = []
        try
            exppost = expect(ψ, "N")
        catch
            exppost = expect(ψ, "Ndn")
        end 
        #corr = correlation_matrix(ψ, "Cdag", "C")

        #@show exppre .- exppost
        @show exppost[ddsite: ddsite + 1]
        #@show corr[ddsite : ddsite + 1, ddsite: ddsite + 1]

        process = Supplywf(ψ)

        workflag = "_repeat$(cnt)"

        run_DPT_many_body(U, L, R,  t_fin; tswitch = tswitch, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  workflag = workflag, ddposition=ddposition,  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "override", vs = vs, process = process, #initlinkdim = initlinkdim
        #n1penalty = α
        stagetype = stagetype,
        TLS = TLS,
        QN = QN
        )

        newα, _ = get_d1QPC(workflag, L, ddsite)

        if abs(newα - α) < 1e-4
            break
        end 

        α = newα
    end 

    


end 



# function DPT_bare()



#     dpt_in = load_JSON( pwd() * "/dptpara.json")
#     U = get(dpt_in, "U", 0.1)
#     L = get(dpt_in, "L", 34)
#     R = get(dpt_in, "R", L)
#     tswitch = Float64(get(dpt_in, "tswitch", 0.0))
#     t_fin = Float64(get(dpt_in, "tfin", 34.0))
#     τ = get(dpt_in, "timestep", 0.25)
#     TEdim = get(dpt_in, "TEdim", 64)
#     biasLR = get(dpt_in, "biasLR", 0.0)
#     mixed = get(dpt_in, "mixed", true)
#     ordering = get(dpt_in, "ordering", "SORTED")
#     avg = get(dpt_in, "avg", false)
#     sweepcnt = get(dpt_in, "sweepcnt", 30)
#     vs = get(dpt_in, "vs", 0.25)


#     workflag = "bare"

#     run_DPT_many_body(U, L, R,  t_fin; tswitch = tswitch, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  workflag = workflag, ddposition="R",  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "disconnectDD", vs = vs)

    


# end 



# function DPT_trend()


#     dpt_in = load_JSON( pwd() * "/dptpara.json")
#     U_start = get(dpt_in, "U", 0.1)
#     L = get(dpt_in, "L", 34)
#     R = get(dpt_in, "R", L)
#     tswitch = Float64(get(dpt_in, "tswitch", 0.0))
#     t_fin = Float64(get(dpt_in, "tfin", 34.0))
#     τ = get(dpt_in, "timestep", 0.25)
#     TEdim = get(dpt_in, "TEdim", 64)
#     biasLR = get(dpt_in, "biasLR", 0.0)
#     mixed = get(dpt_in, "mixed", true)
#     ordering = get(dpt_in, "ordering", "SORTED")
#     avg = get(dpt_in, "avg", false)
#     sweepcnt = get(dpt_in, "sweepcnt", 30)
#     vs = get(dpt_in, "vs", 0.25)
#     fitmode = get(dpt_in, "fitmode", "sqrt")

#     stepsize = 0.1

#     # we first check if there are already existing values
#     existUs = checkunique(glob("work_*", pwd()), "U", x -> parse(Float64, x))

#     if length(existUs) > 2

#         if fitmode == "sqrt"
#             d1s = [ sqrt(get_d1QPC( "_U$(U)", L)[1]) for U in existUs[2:3]]
#         elseif fitmode == "linear"
#             d1s = [ get_d1QPC( "_U$(U)", L)[1] for U in existUs[2:3]]
#         end 

#         U = existUs[1]

#         @warn "RESTART from U = $(U)"

#     else
#         d1s = []
#         Us = [ U_start - stepsize, U_start]
#         # establish initial guess

#         for U in Us

#             workflag = "_U" * string(trunc(U, sigdigits=5)) 
#             run_DPT_many_body(U, L, R,  t_fin; tswitch = tswitch, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  workflag = workflag, ddposition="R",  avg=avg, TEdim = TEdim, sweepcnt = sweepcnt, mode = "disconnectDD", vs = vs)

#             nd1, _ = get_d1QPC(workflag, L)

#             if fitmode == "sqrt"
#                 append!(d1s, sqrt(nd1))
#             elseif fitmode == "linear"
#                 append!(d1s, nd1)
#             end 

#         end 

#         U = U_start - stepsize * 2
#     end 
        
    
#     while U > 0.01

#         U = trunc(U, sigdigits=5)

#         if fitmode == "sqrt"
#             α = (max(1/sqrt(2), d1s[1] - (d1s[2] - d1s[1])))^2
#         elseif fitmode == "linear"
#             α = max(1/2, d1s[1] - (d1s[2] - d1s[1]))
#         else
#             error("unknown fitmode")
#         end 

#         @show d1s, α
#         workflag = "zero_U$(U)"

#         ψ = run_DPT_many_body(U, L, R,  0.0; tswitch = 0.0, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  ddposition="R",  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "manualempty", n1init = α, vs = vs, ordering = ordering, workflag = workflag, initdd = "EMPTY")
        
#         ψ = two_site_rotate(ψ, α, 0)
#         expval = expect(ψ, "N")
#         corr = correlation_matrix(ψ, "Cdag", "C")

#         @show expval[end - 1: end]
#         @show corr[end - 1: end, end - 1: end]

#         process = Supplywf(ψ)

#         workflag = "_U$(U)"

#         run_DPT_many_body(U, L, R,  t_fin; tswitch = tswitch, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  workflag = workflag, ddposition="R",  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "override", vs = vs, process = process)

#         nd1, _ = get_d1QPC(workflag, L)

#         if fitmode == "sqrt"
#             d1s = [sqrt(nd1), d1s[1]]
#         elseif fitmode == "linear"
#             d1s = [nd1, d1s[1]]
#         end 

#         U -= stepsize
#     end 

# end 
