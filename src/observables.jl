

"""Calculates transition charge density between two wf"""
function cal_TCD(ψ1, ψ2; start_off = 0, end_off = 0
    #temp=true
    )

  # ψ1 is left, ψ2 is right

  #L = length(ψ1)
  # systype? fix if necessary
  #tcd = zeros( ComplexF64, (L) )
  # get site indices

  #s1 = siteinds(ψ1)

  # for j = 1:L
  #   replaceind!(ψ1[j], s1[j], s2[j])
  # end 
    operator = "N"
#   W = MPO( s2,  operator)

#   tcd = inner(ψ1', W, ψ2
    tcd = real(inner_product(ψ1, ψ2, operator))[ 1 + start_off : end - end_off]
    #@show tcd
  # explicitly calculate inner product of post-operated states

#   for i in 1:L
#     #tcd[ i] = inner( ψ1', op( operator, s2, i), ψ2)

#     Op = op( operator, s2, i)
#     Nψ2 = apply( Op, ψ2)

#     #tcd[i] = (lefts[i] * ψdag[i] * cur * rights[i])[]
#     tcd[i] = inner( ψ1', Nψ2)

#   end 

  #print(tcd)

  return tcd
end 


"""calculate reduced density matrix of one site"""
function RDM(ψ::MPS, b::Int)

    s = siteinds(ψ)

    orthogonalize!(ψ, b)

    # ρ = prime(ψ[b], s[b]) * dag(ψ[b])
    # w = eigen(ρ).spec.eigs


    # SvNRDM = 0.0
    # for val in w

    #     p= val
    #     SvNRDM -= p * log(p)
    # end 

    # @show SvNRDM


    U, S, V = svd(ψ[b], (s[b], ))

    Sb = 0.0
    for n in 1:dim(S, 1)
        p = S[n,n]^2
        Sb -= p * log(p)   
    end

    return Sb


end 

""" calculates von neumann and other higher order renyi entropies"""
function entropies(psi::MPS, b::Int, max_order::Int)

    s = siteinds(psi)  
    orthogonalize!(psi, b)

    if b==1
        U,S,V = svd(psi[b], (s[b],))
    else
        U,S,V = svd(psi[b], (linkind(psi, b-1), s[b]))
    end

    SvN = 0.0
    Renyi = [0.0 for _ in 2:max_order]
    for n in 1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)

        for (i, order) in enumerate(2:max_order)
            Renyi[i] += p^order
        end     
        
    end

    for (i, order) in enumerate(2:max_order)
        Renyi[i] = 1/(1 - order) * log(Renyi[i])
    end  
    
    # return as uniform array
    return [SvN, Renyi...]
end



function scan_ee(ψ::MPS, max_order::Int)

    L = length(ψ)
    soi = []

    for s in 1:L 
        append!(soi, s)
    end 

    return vectomat([ entropies(ψ, i, max_order) for i in soi]), [ maximum(size(ψ[ j])) for j in soi] 
  
end 

function scan_RDM(ψ::MPS)

    return [RDM(ψ, b) for b in 1:length(ψ)]

end 




function dyna_SRDM(; ψ=nothing, kwargs...)

    workdir = getworkdir()

    if isnothing(ψ)
        T = []
        SRDMs = []

        for file in get_dyna_files()

            ψ = load_ψ(file)
            t = get_time(file)
            append!(T, t)

            println("Calculating SRDM, $t")

            SRDM = scan_RDM(ψ)
            @show append!(SRDMs, [SRDM])

        end 

        writedlm(workdir* "times", T)
        writedlm(workdir* "SRDM", SRDMs)
    else

        SRDM = scan_RDM(ψ)
        open( workdir * "SRDM", "a") do io
            writedlm(io, [SRDM])
        end 

    end 

end 



function dyna_EE(; max_order=3, ψ=nothing, kwargs...)


    #sites = []
    workdir = getworkdir()

    if isnothing(ψ)
        T = []
        bonds = []
        SvN = []
        SRenyi = [[] for _ in 2:max_order]

        for file in get_dyna_files()

            ψ = load_ψ(file)
            t = get_time(file)

            println("Calculating EE, $t")
            ee, bond = scan_ee(ψ, max_order)
            
            print(ee)
            vN = ee[:, 1]

            for (i, order) in enumerate(2:max_order)
                Renyi = ee[:, order]
                append!(SRenyi[i], [Renyi])
            end 
            
            append!(T, t)
            append!(SvN, [vN])
            append!(bonds, [bond])
            #append!(sites, [site])

        end 

        writedlm(workdir * "times", T)
        writedlm(workdir * "SvN", SvN)
        writedlm(workdir * "bonds", bonds)

        for (i, order) in enumerate(2:max_order)
            writedlm(workdir * "SRenyi" * string(order), SRenyi[i])
        end 
        #writedlm(workdir * "sites", sites)

    else
        ee, bond = scan_ee(ψ, max_order)

        open(workdir * "SvN", "a") do io
            writedlm(io, [ee[:,1]])
        end 

        for order in 2:max_order
            open(workdir * "SRenyi" * string(order), "a") do io
                writedlm(io, [ee[:,order]])
            end 
        end 

        open(workdir * "bonds", "a") do io
            writedlm(io, [bond])
        end 

    end 

end 

function static_occ(;static_str="temp_plasmon")

    @info "calculating static occ"
    occ = []


    for file in get_static_files(static_str)

        ψ = load_ψ(file, tag= "psi")
        append!(occ, [expect(ψ, "N")])

    end 

    open(getworkdir() * "staticocc", "w") do io
        writedlm( io, occ)
    end 

end 

function static_tcd(;padding=false, start_off=0, end_off=0, static_str="temp_plasmon")

    tcd = []

    files = get_static_files(static_str)

    gs = load_ψ( TEMP_tag * static_str * "1",  tag="psi")

    for file in files

        ψ = load_ψ(file, tag= "psi")

        res = cal_TCD(ψ, gs; start_off=start_off, end_off=end_off)

        if padding
            res = vcat([0], res, [0])
            res = (res[1:end-1] + res[2:end])/2
        end 

        append!(tcd, [res])

    end 

    tcd = vectomat(tcd)
    open(getworkdir() * "TCD", "w") do io
        writedlm( io, tcd)
    end 

    return tcd

end 


function dyna_tcd(; gs=get_tcd_gs(), ψ=nothing,  kwargs...)

    workdir = getworkdir()

    sys = get(kwargs, :sys, set_QE_two())

    if typeof(sys) == QE_two
        start_off = end_off = 2
    else
        start_off = end_off = 0
    end 

    if isnothing(ψ)
        T = []
        TCD = []

        for file in get_dyna_files()

            ψ = load_ψ(file)
            t = get_time(file)
            append!(T, t)

            println("Calculating TCD, $t")

            tcd = cal_TCD(ψ, gs; start_off=start_off, end_off=end_off)
            @show append!(TCD, [tcd])

        end 

        writedlm(workdir* "times", T)
        writedlm(workdir* "TCDdyna", TCD)
    else

        tcd = cal_TCD(ψ, gs; start_off=start_off, end_off=end_off)
        open( workdir * "TCDdyna", "a") do io
            writedlm(io, [tcd])
        end 

    end 


end 




function dyna_occ(; sys=set_Chain(), ψ=nothing, kwargs...)

    workdir = getworkdir()

    if isnothing(ψ)
        if systype(sys) == "Electron"
            occup = []
            occdn = []
        else
            occs =[]
        end 

        T = []

        for file in get_dyna_files()

            ψ = load_ψ(file)
            t = get_time(file)
            append!(T, t)

            println("Calculating occ, $t")

            if systype(sys)== "Electron"
                append!(occup, [expect(ψ, "Nup")])
                append!(occdn, [expect(ψ, "Ndn")])
        
            else
                @show occ = expect(ψ, "N")
                append!(occs, [occ])
            end 

        end 

        if systype(sys) == "Electron"
            writedlm(workdir * "occup", occup)   
            writedlm(workdir * "occdn", occdn)   
    
        else
            writedlm(workdir * "occ", occs)   
        end 

        writedlm(workdir* "times", T)

    else

        if systype(sys) == "Fermion"
            open(workdir * "occ", "a") do io
                writedlm( io, [expect(ψ, "N")])
            end 

        else
            open(workdir * "occup", "a") do io
                writedlm( io, [expect(ψ, "Nup")])
            end 

            open(workdir * "occdn", "a") do io
                writedlm( io, [expect(ψ, "Ndn")])
            end 

        end 

    end 

end 

function dyna_coherence(; ψ=nothing, sys=set_DPT(), kwargs...)

    function work(ψ, sys)
        
        if systype(sys) == "Fermion" 
            op1, op2 = "Cdag", "C"
        else
            op1, op2 = "Cdagup", "Cup"
        end 
        
        coh = correlation_matrix(ψ, op1, op2; sites= get_systotal(sys) -1 : get_systotal(sys))
        return coh[1, 2]
    end 

    workdir = getworkdir()

    if isnothing(ψ)
        T = []
        cohs = []

        for file in get_dyna_files()

            ψ = load_ψ(file)
            t = get_time(file)
            append!(T, t)

            println("Calculating coherence, $t")

            coh = work(ψ, sys)
            @show append!(cohs, coh)

        end 

        writedlm(workdir* "times", T)
        writedlm(workdir* "coherence", cohs)
    else

        coh = work(ψ, sys)
        open( workdir * "coherence", "a") do io
            writedlm(io, coh)
        end 

    end 


end 

function dyna_dptcurrent(; ψ=nothing, sys=set_DPT(), kwargs...)

    function work(ψ, sys)
        
        if systype(sys) == "Fermion" 
            op1, op2 = "Cdag", "C"
        else
            op1, op2 = "Cdagup", "Cup"
        end 
        
        corr = correlation_matrix(ψ, op1, op2; sites=L_end(sys):R_begin(sys))
        return 2 * imag(corr[1, end])
    end 

    workdir = getworkdir()

    if isnothing(ψ)
        T = []
        currentLR = []

        for file in get_dyna_files()

            ψ = load_ψ(file)
            t = get_time(file)
            append!(T, t)

            println("Calculating DPT current, $t")

            current = work(ψ, sys)
            @show append!(currentLR, current)

        end 

        writedlm(workdir* "times", T)
        writedlm(workdir* "currentLR", currentLR)
    else

        current = work(ψ, sys)
        open( workdir * "currentLR", "a") do io
            writedlm(io, current)
        end 

    end 


end 

function dyna_dptcurrent_mix(; ψ=nothing, sys=set_DPT_mixed(), kwargs...)

    function work(ψ, ULR, sys)
        # except the DD sites
        corr = correlation_matrix(ψ, "Cdag", "C"; sites=union(L_begin(sys):L_end(sys), R_begin(sys):R_end(sys)))
        
        @assert size(ULR) == size(corr)
        return 2 * imag( sum(ULR .* corr))
    end 


    workdir = getworkdir()

    ks = readdlm( workdir * "ks")
    LR = readdlm( workdir * "LR")

    UL, UR = Uk(1, ks, LR)
    #ULR = UL .* UR'
    ULR = UR .* UL'

    if isnothing(ψ)
        T = []
        currentLR = []

        #@show diag(ULR)

        for file in get_dyna_files()

            ψ = load_ψ(file)
            t = get_time(file)
            append!(T, t)

            println("Calculating DPT mix current, $t")

            current = work(ψ, ULR, sys)
            append!(currentLR, current)

        end 

        writedlm(workdir* "times", T)
        writedlm(workdir* "currentLR", currentLR)
    else

        current = work(ψ, ULR, sys)
        open( workdir * "currentLR", "a") do io
            writedlm(io, current)
        end 
        
    end 


end 

function dyna_SDcurrent(; ψ=nothing, sys :: SD_array=set_SD(), t=nothing, kwargs...)
    function work(ψ, sys)
        
        if  sys.systype == "Fermion" 
            ops = [["Cdag", "C"]]
        else
            ops = [["Cdagup", "Cup"]]
        end 
        
        corrs = corr_work(ψ, ops, t)[1]
        current = []

        source = sys.source
        drain = sys.drain
        arr = sys.array

        offset = get_systotal(source) + get_systotal(arr)

        if typeof(source) == Reservoir_spatial

            for (coupling, site) in source.ext_contact
                contact = source.contact
                @show contact, site
                currentval = -2 * coupling * imag(corrs[ contact, site])
                append!(current, currentval)
            end 

            for (coupling, site) in drain.ext_contact
                contact = drain.contact + offset
                @show contact, site
                currentval = -2 * coupling * imag(corrs[ site, contact])
                append!(current, currentval)

            end 
            # 

        else
            LR = vcat(source.LR, drain.LR)
            ks = vcat(source.ks, drain.ks)

            sourceinds = filter( x -> LR[x] > 0, 1:length(LR))
            draininds = filter( x -> LR[x] < 0, 1:length(LR))

            sourceadjinds = [ x > get_systotal(source) ? x + get_systotal(arr) : x for x in  sourceinds]
            drainadjinds = [ x > get_systotal(source) ? x + get_systotal(arr) : x for x in  draininds]

            Usource = [ Ujk(source, 1, k ) for k in ks[sourceinds]]
            Udrain = [ Ujk(drain, 1, k) for k in ks[draininds]]

            @show sourceinds, draininds, sourceadjinds, drainadjinds

            # the contacts has both the source and drain
            for (coupling, site) in source.ext_contacts[1]
                
                currentval = sum(Usource .* corrs[ sourceadjinds, site])
                currentval = - 2 * coupling * imag(currentval)
                append!(current, currentval)
            end 

            for (coupling, site) in source.ext_contacts[2]

                currentval = sum(Udrain .* corrs[ site, drainadjinds])
                currentval = - 2 * coupling * imag(currentval)
                append!(current, currentval)
            end 

        end 

        # check 
        return current
    end 

    workdir = getworkdir()

    if isnothing(ψ)
        T = []
        SDcurrents = []

        for file in get_dyna_files()

            ψ = load_ψ(file)
            t = get_time(file)
            append!(T, t)

            println("Calculating DPT current, $t")

            current = work(ψ, sys)
            @show append!(SDcurrents, [current])

        end 

        writedlm(workdir* "times", T)
        writedlm(workdir* "currentSD", SDcurrents)
    else

        current = work(ψ, sys)
        open( workdir * "currentSD", "a") do io
            writedlm(io, [current])
        end 

    end 


end 



function corr_work(ψ ::MPS,  ops:: Vector{Vector{String}}, t)

    workdir = getworkdir()
    corrs = []

    for (op1, op2) in ops
        corr = correlation_matrix(ψ, op1, op2)
        #@show corr
        outfile = workdir * "corr" * op1 * op2 * ".h5"
        h5open(outfile, isfile( outfile) ? "r+" : "w") do io

            if !haskey(io, string(t))
                write(io, string(t), corr)
            else
                @warn "duplicate key found, not write"
            end 
            
        end 
        
        append!(corrs, [corr])
    end 

    return corrs

end 


function dyna_corr(; ψ=nothing, sys=set_Chain(), t=nothing, kwargs...) 

    workdir = getworkdir()
    ops = ops_determiner(sys)

    if isnothing(ψ)

        T = []
        
        for file in get_dyna_files()

            ψ = load_ψ(file)
            t = get_time(file)
            println("Calculating corr, $t")
            _ = corr_work(ψ, ops, t)

            append!(T, t)
        end 

        writedlm(workdir* "times", T)

    else
        
        _ = corr_work(ψ, ops, t)
    end 

    return nothing

end 




function dyna_LSRcurrent()

    T = []
    current = []

    workdir = getworkdir()

    for file in get_dyna_files()

        ψ = load_ψ(file)
        t = get_time(file)
        append!(T, t)

        println("Calculating LSR current, $t")

        L = div(length(ψ), 2) 
        println(L)
        corr = correlation_matrix(ψ, "Cdag", "C")
        append!(current, - 2 * (imag(corr[L, L + 1]) - imag(corr[1, 2])))

    end 

    writedlm(workdir* "times", T)
    writedlm(workdir* "current", current)


end 

function dyna_pϕ()

    pϕ = []
    T = []

    workdir = getworkdir()

    for file in get_dyna_files()

        ψ = load_ψ(file)
        t = get_time(file)
        append!(T, t)

        println("Calculating partial contraction, $t")

        #state = [ partial_contract(ψ, [i, i+1]) for i in 1:div(length(ψ), 2) ]
        state = partial_contract(ψ, [1, 2])

        normalize!(state)
        c1, c2 = state[1, 2], state[2, 1]
        ϕ1, ϕ2 = angle(c1), angle(c2)
        p1 = abs2(c1) / (abs2(c1) + abs2(c2))
        ϕ = ((ϕ2 - ϕ1)/pi )% 2

        append!(pϕ, [[p1, ϕ]])

    end 


    writedlm(workdir* "times", T)
    writedlm(workdir* "coefficients", pϕ)

end 