

function QE_static(key, QEen, output, product; TEdim=64, τ=1.0, dp=1.0, staticex= 1, QEmul=1.0, center_parameter = EMPTY_CENTER, kwargs...)

    @info "Begin. QEen = $QEen"
    QEen = get_QEen(QEen, key, output, TEdim, QEmul, product; kwargs...)
    @info "After getting QEen, QEen = $QEen"

    ψ = !product ? load_ψ(output) : gen_state(sys)

    sys = QE_determiner(key; QEen=QEen, dp=dp, center_parameter = center_parameter, kwargs...)
    static = set_Static(; ex=staticex, sweepdim=TEdim, kwargs...)
    run_static_simulation(sys, static, ψ)

end 

function QE_dyna_regular(key, QEen, output, product; TEdim=64, τ=1.0, dp=1.0,  QEmul=1.0, start=τ ,fin=200.0, center_parameter = EMPTY_CENTER, save_every=false, adiabatic=0.0, kwargs...)

    @info "Begin. QEen = $QEen"
    QEen = get_QEen(QEen, key, output, TEdim, QEmul, product; kwargs...)
    @info "After getting QEen, QEen = $QEen"
    

    obs = [dyna_EE, dyna_occ, dyna_corr]
    #start = get(kwargs, :start, τ)
    init_key = start == τ ? output : LASTSTSTR
    ψ = !product ? load_ψ(init_key) : gen_state(sys)

    @assert adiabatic <= tswitch

    if adiabatic != 0 && start < adiabatic + τ

        for t in start:τ:adiabatic

            @show temp_dp = dp * t / adiabatic

            temp_sys = QE_determiner(key; QEen=QEen, dp=temp_dp, center_parameter = center_parameter, kwargs...)
            # we only run it for 1 step
            temp_dynamic = set_Dynamic(; TEdim=TEdim, τ=τ, start=t, fin=t, kwargs..., )
            ψ = run_dynamic_simulation(temp_sys, temp_dynamic, ψ; message="QEdyna_temp", save_every=save_every, obs=obs)

        end 
        
        start = adiabatic + τ
    end 

    sys = QE_determiner(key; QEen=QEen, dp=dp, center_parameter = center_parameter, kwargs...)
    dynamic = set_Dynamic(; TEdim=TEdim, τ=τ, start=start, fin=fin, kwargs...)
    run_dynamic_simulation(sys, dynamic, ψ; message="QEdyna", save_every=save_every, obs=obs)

    
end 


function QE_confine(key, QEen, output, product; TEdim=64, τ=1.0, dp=1.0,  QEmul=1.0, start=τ ,fin=200.0, center_parameter = EMPTY_CENTER, save_every=false, tswitch=0.0, confine_parameters = EMPTY_CONFINES, kwargs...)

    @info "Begin. QEen = $QEen"
    QEen = get_QEen(QEen, key, output, TEdim, QEmul, product; confine_parameters=confine_parameters, kwargs...)
    @info "After getting QEen, QEen = $QEen"
    

    obs = [dyna_EE, dyna_occ, dyna_corr]
    #start = get(kwargs, :start, τ)
    init_key = start == τ ? output : LASTSTSTR
    ψ = !product ? load_ψ(init_key) : gen_state(sys)

    if start < tswitch
        @info "Stage 1"
        Stage1 = QE_determiner(key; QEen=QEen, dp=dp, center_parameter = center_parameter, confine_parameters=confine_parameters, kwargs...)
        dynamic = set_Dynamic(; TEdim=TEdim, τ=τ, start=start, fin=tswitch, kwargs...)
        ψ = run_dynamic_simulation(Stage1, dynamic, ψ; message="QEStage1", save_every=save_every, obs=obs)

    else
        @warn "tswitch <= start, is this the expected behavior?"
        tswitch = start
    end 

    #Stage 2 we release the potential

    @info "Stage 2"
    Stage2= QE_determiner(key; QEen=QEen, dp=dp, center_parameter = center_parameter, confine_parameters=EMPTY_CONFINES, kwargs...)
    dynamic = set_Dynamic(; TEdim=TEdim, τ=τ, start=tswitch + τ, fin=fin, kwargs...)
    run_dynamic_simulation(Stage2, dynamic, ψ; message="QEStage2", save_every=save_every, obs=obs)


    
end 


function QE_wrapper(key)

    product = false
    qe_in = load_JSON( pwd() * "/qepara.json")
    QEen = get(qe_in, "QEen", 0.0)
    L = get(qe_in, "L", 20)
    N = get(qe_in, "N", div(L, 2))
    QEmul = get(qe_in, "QEmul", 1.0)
    dp = get(qe_in, "dp", 1.0)

    #τ = get(qe_two_in, "timestep", 0.125)
    TEdim = get(qe_in, "TEdim", 64)
    center_parameter = get(qe_in, "center_parameter", EMPTY_CENTER)
    confine_parameters = get(qe_in, "confine_parameters", EMPTY_CONFINES)

    τ = get(qe_in, "timestep", 1.0)
    start = get(qe_in, "start", τ)
    fin = get(qe_in, "fin", 200.0)
    tswitch = get(qe_in, "tswitch", 50.0)
    #adiabatic = get(qe_in, "adiabatic", 0.0)

    inits = get(qe_in, "inits", "1")

    
    #run_QE_two(QEen, L, N, product; staticex= 0, dp=1.0, QEmul=QEmul, TEdim=TEdim)
    # QE_dyna_regular("QE_HOM", QEen, "initialqeparallelstate", product; QEmul=QEmul, TEdim=TEdim, L=L, N=N, τ=τ,  start = start, fin=fin, center_parameter=center_parameter, dp=dp, inits=inits, adiabatic=adiabatic)
    QE_confine(key, QEen, "initialqeparallelstate", product; QEmul=QEmul, TEdim=TEdim, L=L, N=N, τ=τ,  start = start, fin=fin, center_parameter=center_parameter, confine_parameters=confine_parameters, dp=dp, inits=inits, tswitch=tswitch)


end 




function QE_gaussian_wrapper()

    qe_gaussian_in = load_JSON( pwd() * "/qegaussian.json")

    full_size = get(qe_gaussian_in, "fullsize", 100)
    full_N = get(qe_gaussian_in, "fullN", div(full_size, 2))
    τ = get(qe_gaussian_in, "timestep", 0.25)
    fin = get(qe_gaussian_in, "fin", 200)
    TEdim = get(qe_gaussian_in, "TEdim", 128)
    center = get(qe_gaussian_in, "center", 6.5)
    sigma = get(qe_gaussian_in, "sigma", 4)
    includegs = get(qe_gaussian_in, "includegs", false)
    padding = get(qe_gaussian_in, "padding", false)
    L = get(qe_gaussian_in, "L", 12)
    conv = get(qe_gaussian_in, "conv", true)
    mode = get(qe_gaussian_in, "mode", "biasedchain")

    static_str = get_static_str(mode)

    save_every=false
    obs = [dyna_EE, dyna_occ, dyna_corr, dyna_tcd]

    # empty static wf, solve for QE
    if isempty(get_static_files(static_str))

        @warn "no static files, solving for static chain eigenstates"
        solve_QE(; para_in = qe_gaussian_in, mode=mode)
    end 
    
    ψ = prepare_wavepacket(; includegs=includegs, center=center, sigma=sigma, L=full_size, padding=padding, conv=conv, mode=mode)

    if mode == "biasedchain"
        sys= set_Chain(; L=full_size, N = full_N )
    elseif mode == "QEtwo"
        sys = set_QE_two(; L=L, N=N)
    end 

    dynamic = set_Dynamic(; τ=τ, TEdim=TEdim, start=τ, fin=fin)
    run_dynamic_simulation(sys, dynamic, ψ; message="QE_gaussian_chain_only", save_every=save_every, obs=obs)

    return nothing

end 

# function QE_SIAM_wrapper()

    
#     prod = false
#     qe_siam_in = load_JSON( pwd() * "/qesiampara.json")

#     QEen = get(qe_siam_in, "QEen", 0.0)
#     legleft = get(qe_siam_in, "legleft", 2)
#     legright = get(qe_siam_in, "legright", legleft)
#     siteseach = get(qe_siam_in, "siteseach", 2)
#     N = get(qe_siam_in, "N", div(siteseach, 2))
#     QEmul = get(qe_siam_in, "QEmul", 1.0)
#     init = get(qe_siam_in, "init", "1")
#     TTN = get(qe_siam_in, "TTN", false)
#     center_parameter = get(qe_siam_in, "center_parameter", EMPTY_CENTER)


#     #τ = get(qe_siam_in, "timestep", 0.125)
#     TEdim = get(qe_siam_in, "TEdim", 64)
#     timestep = get(qe_siam_in, "timestep", 1.0)

#     run_QE_SIAM("QE_SIAM", QEen, "initialqesiamstate", prod; siteseach=siteseach, N=N,  TTN=TTN, init= init, legleft=legleft, legright=legright, τ=timestep, QEmul=QEmul, TEdim = TEdim, center_parameter=center_parameter)

#     dyna_occ()
#     dyna_EE()
#     #dyna_dptcurrent()


    

# end 