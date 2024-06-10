

function run_QE(key, QEen, output, product; TEdim=64, τ=1.0, dp=1.0, staticex= 0, QEmul=1.0, start=τ ,fin=200, center_parameter = EMPTY_CENTER, save_every=false, adiabatic=0.0, kwargs...)

    try
        QEen = load_plsmon(output) * QEmul
        println("QE file found, QE energy = ", QEen)
    catch
        println("QE file not found,beginning next stage, QE energy = ", QEen)
    end


    # if no QEen, we need to perform further calculations on the initial state
    # the basic logic is that this is a QE calculation, QEen =0 makes no sense
    if QEen == 0 || (!product && !check_ψ(output))
        

        println("Calculating initial state")
        ex = (QEen == 0) ? 2 : 1
        # we set up the decoupled sys from the QE
        decoupled = QE_determiner(key; QEen=0.0, dp=0.0, center_parameter = EMPTY_CENTER, kwargs...)

        # get plasmon energy
        static = set_Static(; ex=ex, output=output, sweepdim=TEdim, kwargs...)

        ϕ = gen_state(decoupled)
        run_static_simulation(decoupled, static, ϕ; message="QEinit")

        QEen = load_plsmon(output) * QEmul

    end 

    
    if staticex == 0

        obs = [dyna_EE, dyna_occ, dyna_corr]
        #start = get(kwargs, :start, τ)
        fin = fin * τ
        init_key = start == τ ? output : LASTSTSTR
        ψ = !product ? load_ψ(init_key) : gen_state(sys)

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


    else
        ψ = !product ? load_ψ(output) : gen_state(sys)

        sys = QE_determiner(key; QEen=QEen, dp=dp, center_parameter = center_parameter, kwargs...)
        static = set_Static(; ex=staticex, sweepdim=TEdim, kwargs...)
        run_static_simulation(sys, static, ψ)
    end 
    

end 



function QE_two_wrapper()

    product = false
    qe_two_in = load_JSON( pwd() * "/qetwopara.json")

    QEen = get(qe_two_in, "QEen", 0.0)
    L = get(qe_two_in, "L", 20)
    N = get(qe_two_in, "N", div(L, 2))
    QEmul = get(qe_two_in, "QEmul", 1.0)

    #τ = get(qe_two_in, "timestep", 0.125)
    TEdim = get(qe_two_in, "TEdim", 64)

    #run_QE_two(QEen, L, N, product; staticex= 0, dp=1.0, QEmul=QEmul, TEdim=TEdim)
    run_QE("QE_two", QEen, "initialqe2state", product; QEmul=QEmul, TEdim=TEdim, L=L, N=N)
    dyna_occ()
    dyna_EE()

end 

function QE_parallel_wrapper()

    product = false
    qe_parallel_in = load_JSON( pwd() * "/qeparallelpara.json")

    QEen = get(qe_parallel_in, "QEen", 0.0)
    L = get(qe_parallel_in, "L", 20)
    N = get(qe_parallel_in, "N", div(L, 2))
    QEmul = get(qe_parallel_in, "QEmul", 1.0)
    dp = get(qe_parallel_in, "dp", 1.0)

    #τ = get(qe_two_in, "timestep", 0.125)
    TEdim = get(qe_parallel_in, "TEdim", 64)
    center_parameter = get(qe_parallel_in, "center_parameter", EMPTY_CENTER)

    τ = get(qe_parallel_in, "timestep", 1.0)
    start = get(qe_parallel_in, "start", τ)
    fin = get(qe_parallel_in, "fin", 200)
    adiabatic = get(qe_parallel_in, "adiabatic", 0.0)

    inits = get(qe_parallel_in, "inits", "1")
    #run_QE_two(QEen, L, N, product; staticex= 0, dp=1.0, QEmul=QEmul, TEdim=TEdim)
    run_QE("QE_HOM", QEen, "initialqeparallelstate", product; QEmul=QEmul, TEdim=TEdim, L=L, N=N, τ=τ,  start = start, fin=fin, center_parameter=center_parameter, dp=dp, inits=inits, adiabatic=adiabatic)
    dyna_occ()
    dyna_EE()

end 


function QE_SIAM_wrapper()

    
    prod = false
    qe_siam_in = load_JSON( pwd() * "/qesiampara.json")

    QEen = get(qe_siam_in, "QEen", 0.0)
    legleft = get(qe_siam_in, "legleft", 2)
    legright = get(qe_siam_in, "legright", legleft)
    siteseach = get(qe_siam_in, "siteseach", 2)
    N = get(qe_siam_in, "N", div(siteseach, 2))
    QEmul = get(qe_siam_in, "QEmul", 1.0)
    init = get(qe_siam_in, "init", "1")
    TTN = get(qe_siam_in, "TTN", false)
    center_parameter = get(qe_siam_in, "center_parameter", EMPTY_CENTER)


    #τ = get(qe_siam_in, "timestep", 0.125)
    TEdim = get(qe_siam_in, "TEdim", 64)
    timestep = get(qe_siam_in, "timestep", 1.0)

    run_QE_SIAM("QE_SIAM", QEen, "initialqesiamstate", prod; siteseach=siteseach, N=N,  TTN=TTN, init= init, legleft=legleft, legright=legright, τ=timestep, QEmul=QEmul, TEdim = TEdim, center_parameter=center_parameter)

    dyna_occ()
    dyna_EE()
    #dyna_dptcurrent()


    

end 