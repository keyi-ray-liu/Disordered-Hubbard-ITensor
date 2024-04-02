function run_QE_two(QEen, dims, N, product; staticex= 0, dp=1.0, kwargs...)

    coulomb = set_Coulombic(;kwargs...)
    chain = set_Chain(;dims=dims, N=N, coulomb=coulomb, kwargs...)

    output = "initialstate"

    # if no QEen, we need to perform further calculations on the initial state
    if QEen == 0 || !product
        
        println("Calculating initial state")
        ex = (QEen == 0) ? 2 : 1
        # we set up the decoupled sys from the QE
        decoupled = set_QE_two(; chain_only=chain, QEen=0.0, dp=0.0, kwargs...)

        # get plasmon energy
        static = set_Static(; ex=ex, output=output, kwargs...)

        ϕ = gen_state(decoupled)
        run_static_simulation(decoupled, static, ϕ)

        QEen = (QEen == 0) ? load_plsmon(output) * QEmul(decoupled) : QEen

    end 

    sys = set_QE_two(; QEen=QEen, chain_only=chain, dp=dp, kwargs...)

    ψ = (!product && staticex == 0) ? load_ψ(output) : gen_state(sys)

    if staticex == 0
        dynamic = set_Dynamic(;kwargs...)
        run_dynamic_simulation(sys, dynamic, ψ)
    else
        static = set_Static(; ex=staticex, kwargs...)
        run_static_simulation(sys, static, ψ)
    end 
    

end 
