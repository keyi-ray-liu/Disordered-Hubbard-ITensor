function run_QE_two(QEen, dims, N, product; staticex= 0, dp=1.0, TEdim=64, τ=1.0, QEmul=1.0, kwargs...)

    coulomb = set_Coulombic(;kwargs...)
    chain = set_Chain(;dims=dims, N=N, coulomb=coulomb, kwargs...)

    output = "initialqe2state"

    # if no QEen, we need to perform further calculations on the initial state
    if QEen == 0 || !product
        
        println("Calculating initial state")
        ex = (QEen == 0) ? 2 : 1
        # we set up the decoupled sys from the QE
        decoupled = set_QE_two(; chain_only=chain, QEen=0.0, dp=0.0, kwargs...)

        # get plasmon energy
        static = set_Static(; ex=ex, output=output, sweepdim=TEdim, kwargs...)

        ϕ = gen_state(decoupled)
        run_static_simulation(decoupled, static, ϕ)

        QEen = (QEen == 0) ? load_plsmon(output) * QEmul : QEen

    end 

    sys = set_QE_two(; QEen=QEen, chain_only=chain, dp=dp, kwargs...)

    ψ = (!product && staticex == 0) ? load_ψ(output) : gen_state(sys)

    if staticex == 0
        dynamic = set_Dynamic(; TEdim=TEdim, τ=τ, start=τ, fin=200*τ, kwargs...)
        run_dynamic_simulation(sys, dynamic, ψ)
    else
        static = set_Static(; ex=staticex, sweepdim=TEdim, kwargs...)
        run_static_simulation(sys, static, ψ)
    end 
    

end 


function run_QE_SIAM(QEen, siteseach, N, product; TTN = false, staticex= 0, TEdim=64, τ=1.0, QEmul=1.0, init="1", kwargs...)

    output = "initialqesiamstate"

    # if no QEen, we need to perform further calculations on the initial state
    if QEen == 0 || !product
        
        println("Calculating initial state")
        ex = (QEen == 0) ? 2 : 1
        # we set up the decoupled sys from the QE
        decoupled = set_QE_SIAM(;  TTN=TTN, init=init, QEen=0.0, dp=0.0, N = N, siteseach=siteseach, kwargs...)

        # get plasmon energy
        static = set_Static(; ex=ex, output=output, sweepdim=TEdim, kwargs...)

        ϕ = gen_state(decoupled)
        run_static_simulation(decoupled, static, ϕ)

        QEen = (QEen == 0) ? load_plsmon(output) * QEmul : QEen

    end 

    sys = set_QE_SIAM(; TTN=TTN, init=init, QEen=QEen, N=N, siteseach=siteseach, kwargs...)

    ψ = (!product && staticex == 0) ? load_ψ(output) : gen_state(sys)

    if staticex == 0
        dynamic = set_Dynamic(; τ=τ, start=τ, fin=200*τ, TEdim=TEdim, kwargs...)
        run_dynamic_simulation(sys, dynamic, ψ)
    else
        static = set_Static(; ex=staticex, kwargs...)
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

    run_QE_two(QEen, L, N, product; staticex= 0, dp=1.0, QEmul=QEmul, TEdim=TEdim)

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


    #τ = get(qe_siam_in, "timestep", 0.125)
    TEdim = get(qe_siam_in, "TEdim", 64)
    timestep = get(qe_siam_in, "timestep", 1.0)

    run_QE_SIAM(QEen, siteseach, N, prod; TTN=TTN, init= init, legleft=legleft, legright=legright, τ=timestep, QEmul=QEmul, TEdim = TEdim)

    dyna_occ()
    dyna_EE()
    #dyna_dptcurrent()


    

end 