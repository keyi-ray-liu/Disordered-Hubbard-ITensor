get_static_str(mode::String) = mode == "QE_two" ?  "QEplasmon" : "plasmon"

function get_QEen(QEen, key, output, TEdim, QEmul, product; kwargs...)
    @info "Begin. QEen = $QEen"

    try
        QEen = load_plasmon(output) * QEmul
        @info "QE file found, QE energy = $QEen"
    catch
        @info "QE file not found, beginning next stage, QE energy = $QEen"
    end

    # if no QEen, we need to perform further calculations on the initial state
    # the basic logic is that this is a QE calculation, QEen =0 makes no sense
    if QEen == 0 || (!product && !check_ψ(output))
        println("Calculating plasmon state")
        ex = (QEen == 0) ? 2 : 1
        # we set up the decoupled sys from the QE
        decoupled = QE_determiner(key; QEen=0.0, dp=0.0, center_parameter = EMPTY_CENTER, kwargs...)
    
        # get plasmon energy
        static = set_Static(; ex=ex, output=output, sweepdim=TEdim, kwargs...)
    
        ϕ = gen_state(decoupled)
        run_static_simulation(decoupled, static, ϕ; message="QEinit")
    
        QEen = load_plasmon(output) * QEmul
    end

    @info "After getting QEen, QEen = $QEen"

    return QEen


end 

function get_QEinit(init_key, key, TEdim; kwargs...)
    

    # if no QEen, we need to perform further calculations on the initial state
    # the basic logic is that this is a QE calculation, QEen =0 makes no sense
    if !check_ψ(init_key)
        @info "calculating separate init state for QE"
        ex = 1
        # we set up the decoupled sys from the QE
        decoupled = QE_determiner(key; QEen=0.0, dp=0.0, center_parameter = EMPTY_CENTER, kwargs...)
    
        # get plasmon energy
        static = set_Static(; ex=ex, output=init_key, sweepdim=TEdim, kwargs...)
    
        ϕ = gen_state(decoupled)
        run_static_simulation(decoupled, static, ϕ; message="QEinit")
    end

    ψ = load_ψ(init_key)
    return ψ


end 




function get_tol_coeff(TCD, L, g, tols; write=false)

    maxvals = [ maximum(abs.(f)) for f in eachrow(TCD)]
    idxs = [findall(x -> x > tol, maxvals) for tol in tols]

    idxs = filter(x->length(eachrow(x)) > 5, idxs)

    Cs = zeros(length(idxs), length(maxvals))

    for (i, idx) in enumerate(idxs)
        F = transpose(TCD[idx, 1:L])
        C = F \ g

        # @show C
        # @show g
        # @show F * C

        #append!(Cs, [C])
        Cs[i, idx] = C

    end 


    if write

        open( getworkdir() * "TCDcoeff", "w") do io
            writedlm(io, Cs)
        end 
        
        open( getworkdir() * "TCDtols", "w") do io
            writedlm(io, tols[1:length(idxs)])
        end 


        open( getworkdir() * "TCDidx", "w") do io
            writedlm(io, idxs)
        end 
    end 

    return Cs, idxs

end 


function wave_coeff(TCD; L=12, center=6.5, sigma=2, conv=false, includegs=true)

    #@show size(TCD[:, 1:12])
    g = map( x-> exp( -((x - center)/sigma) ^ 2), 1:L)

    #@show start_idx = 2 - Int(includegs)

    if conv
        g = g .* TCD[2, 1:L]
    end 

    tols = 1e-4:2e-3:1e-1
    _ = get_tol_coeff(TCD, L, g, tols; write=true)

    Cs, idxs = get_tol_coeff(TCD, L, g, [0.0261])

    return Cs[idxs[1]], idxs[1]
end 

function prepare_wavepacket(; includegs=false, center=6.5, sigma=2, L=12, padding=false, conv=false, mode="QE_two")

    static_str = get_static_str(mode)

    #static_occ(;static_str=static_str)

    if !isfile( getworkdir() * "initwavepacket.h5")

        @warn "no initial wavepacket, calculating"
        if isfile( getworkdir() * "TCD")
            TCD = readdlm( getworkdir() * "TCD")
        else
            @warn "no TCD, calculating"

            if mode == "QE_two"
                start_off = end_off = 2
            else
                start_off = end_off = 0
            end 

            TCD = static_tcd(; padding=padding, start_off=start_off, end_off=end_off, static_str=static_str)
        end 

        coeff, idx = wave_coeff(TCD; center=center, sigma=sigma, L=L, conv=conv, includegs=includegs)

        @show coeff, idx

        wf = [ load_ψ(f; tag="psi") for f in get_static_files(static_str)]
        
        
        #@show start_idx = 2 - Int(includegs)
        wf_to_sum = coeff .* wf[idx]

        ψ = add(wf_to_sum...; maxdim=256)

        normalize!(ψ)

        open(getworkdir() * "combTCD", "w") do io
            writedlm(io, cal_TCD(ψ, wf[1]))
        end 

        h5open( getworkdir() * "initwavepacket.h5", "w") do io
            write(io, "psi", ψ)
        end 

    end 

    @info "loading initial wavepacket"
    ψ = load_ψ("initwavepacket"; tag="psi")
    return ψ

end 


function solve_QE(; para_in = nothing, mode="biased_chain", output=get_static_str(mode), kwargs...)

    if !(typeof(para_in) <: Dict)

        @warn "no initial chain para dict, loading biasedchain"
        para_in = load_JSON( pwd() * "/qegaussian.json")
    end 
    
    @info "QE solver"
    @show para_in
    
    L = get(para_in, "L", 20)
    N = get(para_in, "N", div(L, 2))

    #τ = get(qe_two_in, "timestep", 0.125)
    dim = get(para_in, "dim", 64)
    ex = get(para_in, "ex", 10)
    sweepcnt = get(para_in, "sweepcnt", 10)

    #static_str = get_static_str(mode)

    #run_chain(L, N, ex; dim=dim)

    if mode == "biasedchain"

        full_size = get(para_in, "fullsize", 100)
        chain_start = get(para_in, "chain_start", 1)

        _ = run_biased_chain(full_size, L, N, ex; dim=dim, sweepcnt=sweepcnt, chain_start = chain_start, output=output, kwargs...)

    elseif mode == "QE_two"

        sys = set_QE_two(; dp=0.0, L=L, N=N)
        static = set_Static(; ex=ex, sweepcnt=sweepcnt, sweepdim=dim, output=output)
        ψ = gen_state(sys)

        @show sys
        _ = run_static_simulation(sys, static, ψ)

    else
        error("Unknown mode!")
    end 

end 


"""Here we attempt to attach Quantum emitters to the system without breaking

We strictly set the following:
sys is the input system WITH QE
ψ will be the input wavefunction

We always generate a NEW site basis and state
"""
function QE_embedding(sys::systems, ψ)

    if typeof(sys) == QE_two

        if length(ψ) + 2 * QESITES != get_systotal(sys)
            error("system size mismatch! ψ length != sys - 2 * QE ")
        end 

        ψn = gen_state(sys)
        ψn.data[1 + QESITES : end - QESITES] = ψ.data
        
    end 

    return [ψn]
end 




function test_embedding_wrapper()

    function save(ψ, name)

        open( "temp/" * name, "w") do io
            writedlm(io, ψ)
        end 
    
        open( "temp/data" * name, "w") do io
            writedlm(io, ψ[1])
        end 
    
    end 

    function compare(ψ1, ψ2)

        for i in 1:3
            print(ψ2.data[i])
        end 

    end 

    ori = set_QE_two(;L=12, N=6, dp=0.0)
    static = set_Static(; ex=1, sweepcnt=5)
    ψ = gen_state(ori)
    ψori = run_static_simulation(ori, static, ψ)

    save(ψori, "native")

    ψ = run_chain(12, 6, 1; sweepcnt=10)[1]

    #ψ = load_ψ("temp_plasmon1"; tag = "psi")
    sys = set_QE_two(;L=12, N=6)

    ψn =  QE_embedding(sys, ψ)

    save(ψn, "embed")
    dynamic = set_Dynamic()

    #compare(ψori[1], ψn[1])

    #run_dynamic_simulation(sys, dynamic, ψn)
end 