

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


function RDM2(ket::MPS, a::Int, b::Int)

    orthogonalize!(ket,a)
    sites = siteinds(ket)
    bra = prime(dag(ket),linkinds(ket))
    
    rho = prime(ket[a],linkinds(ket,a-1)) * prime(bra[a], sites[a])
    for j in a+1:b-1
        rho *= ket[j]*bra[j]
    end
    rho *= prime(ket[b],linkinds(ket,b)) * prime(bra[b], sites[b])  

    return rho

end 





"""calculate reduced density matrix of one site"""
function RDM1(ψ::MPS, b::Int)

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


""" calculates explicitly the SVD"""
function SVD(psi::MPS, b::Int)

    s = siteinds(psi)  
    orthogonalize!(psi, b)

    if b==1
        _,S,_ = svd(psi[b], (s[b],))
    else
        _,S,_ = svd(psi[b], (linkind(psi, b-1), s[b]))
    end

    Sval = [ S[n, n] for n in 1:dim(S, 1)]
    Sval = reverse(sort(Sval))
    # return as uniform array
    #@show size(Sval)

    # S = diag(dense(S))
    # Sval = [ val for val in S]
    return Sval
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

    return [RDM1(ψ, b) for b in 1:length(ψ)]

end 



function dyna_SVD(; ψ=nothing, workflag = "", curtime= nothing, kwargs...)

    if isnothing(ψ)
        error("SVD currently only supports in time cal")
    end 

    workdir = getworkdir(workflag)

    S = [SVD(ψ, i) for i in eachindex(ψ)]
    
    #S = resize!.(S, maximum(length, S))
    
    newS = zeros( size(S)[1], maximum(length, S))
    for (i, s) in enumerate(S)
        newS[i, 1: length(s)] = s
    end 


    outfile = workdir * "SVD.h5"

    if curtime % 5 == 0
        h5open( outfile, isfile( outfile) ? "r+" : "w") do io

            if !haskey(io, string(curtime))
                write(io, string(curtime), newS)
            else
                @warn "duplicate key found, not write"
            end 
            
        end 
    end 

end 



function dyna_SRDM(; ψ=nothing, workflag = "", kwargs...)

    workdir = getworkdir(workflag)

    if isnothing(ψ)
        T = []
        SRDMs = []

        for file in get_dyna_files()

            ψ = load_ψ(file, workflag)
            curtime = get_time(file)
            append!(T, curtime)

            println("Calculating SRDM, $curtime")

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

function dyna_product_state_overlap(; ψ=nothing, sys=nothing, workflag = "", ordering = "SORTED", kwargs...)

    workdir = getworkdir(workflag)
    

    if isnothing(ψ)
        error("No supported")
    else

        if typeof(sys.source) == Reservoir_momentum
            sites = siteinds(ψ)
            ϕ = fermilevel(sys, sites, ordering)
            overlap = abs2(inner(ϕ, ψ))

            open( workdir * "productstateoverlap", "a") do io
                writedlm(io, overlap)
            end 
        end 

    end 
    
end 


function dyna_EE(; max_order=1, ψ=nothing, workflag = "", kwargs...)


    #sites = []
    workdir = getworkdir(workflag)

    if isnothing(ψ)
        T = []
        bonds = []
        SvN = []
        SRenyi = [[] for _ in 2:max_order]

        for file in get_dyna_files()

            ψ = load_ψ(file, workflag)
            curtime = get_time(file)

            println("Calculating EE, $curtime")
            ee, bond = scan_ee(ψ, max_order)
            
            print(ee)
            vN = ee[:, 1]

            for (i, order) in enumerate(2:max_order)
                Renyi = ee[:, order]
                append!(SRenyi[i], [Renyi])
            end 
            
            append!(T, curtime)
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
            writedlm(io, [round.(ee[:,1], sigdigits=6)])
        end 

        for order in 2:max_order
            open(workdir * "SRenyi" * string(order), "a") do io
                writedlm(io, [ee[:,order]])
            end 
        end 

        open(workdir * "bonds", "a") do io
            writedlm(io, [round.(bond, sigdigits=6)])
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

    gs = load_ψ( TEMP_tag * static_str * "1";  tag="psi")

    for file in files

        ψ = load_ψ(file; tag= "psi")

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

    sys = get(kwargs, :sys, QE_two())

    if typeof(sys) == QE_two
        start_off = end_off = 2
    else
        start_off = end_off = 0
    end 

    if isnothing(ψ)
        T = []
        TCD = []

        for file in get_dyna_files()

            ψ = load_ψ(file, workflag)
            curtime = get_time(file)
            append!(T, curtime)

            println("Calculating TCD, $curtime")

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




function dyna_occ(; sys=Chain(), ψ=nothing, workflag = "", kwargs...)

    workdir = getworkdir(workflag)

    if isnothing(ψ)
        if systype(sys) == "Electron"
            occup = []
            occdn = []
        else
            occs =[]
        end 

        T = []

        for file in get_dyna_files()

            ψ = load_ψ(file, workflag)
            curtime = get_time(file)
            append!(T, curtime)

            println("Calculating occ, $curtime")

            if systype(sys)== "Electron"
                append!(occup, [expect(ψ, "Nup")])
                append!(occdn, [expect(ψ, "Ndn")])
        
            else
                @show occ = expect(ψ, "N")
                append!(occs, [occ])
            end 

        end 

        if systype(sys) == "Electron"
            writedlm(workdir * "occup", round.(occup, digits=6))   
            writedlm(workdir * "occdn", round.(occdn, digits=6))   
    
        else
            writedlm(workdir * "occ", round.(occs, digits=6))   
        end 

        writedlm(workdir* "times", T)

    else

        if systype(sys) == "Fermion"
            open(workdir * "occ", "a") do io
                writedlm( io, [round.(expect(ψ, "N"), sigdigits=6)])
            end 

        else
            open(workdir * "occup", "a") do io
                writedlm( io, [round.(expect(ψ, "Nup"), sigdigits=6)])
            end 

            open(workdir * "occdn", "a") do io
                writedlm( io, [round.(expect(ψ, "Ndn"), sigdigits=6)])
            end 

        end 

    end 

end 

function dyna_coherence(; ψ=nothing, sys=DPT(), workflag = "", kwargs...)

    function work(ψ, sys)
        
        if systype(sys) == "Fermion" 
            op1, op2 = "Cdag", "C"
        else
            op1, op2 = "Cdagup", "Cup"
        end 

        ddsite = dd_lower(sys)
        
        coh = correlation_matrix(ψ, op1, op2; sites= ddsite : ddsite + 1)

        return coh[1, 2]
    end 

    workdir = getworkdir(workflag)

    if isnothing(ψ)
        T = []
        cohs = []

        for file in get_dyna_files()

            ψ = load_ψ(file, workflag)
            curtime = get_time(file)
            append!(T, curtime)

            println("Calculating coherence, $curtime")

            coh = work(ψ, sys)
            @show append!(cohs, coh)

        end 

        writedlm(workdir* "times", T)
        writedlm(workdir* "coherence", cohs)
    else

        coh = work(ψ, sys)
        open( workdir * "coherenceRE", "a") do io
            writedlm(io, real.(coh))
        end 

        open( workdir * "coherenceIM", "a") do io
            writedlm(io, imag.(coh))
        end 

    end 


end 

function dyna_dptcurrent(; ψ=nothing, sys=DPT(), workflag = "", kwargs...)

    function work(ψ, sys)
        
        if systype(sys) == "Fermion" 
            op1, op2 = "Cdag", "C"
        else
            op1, op2 = "Cdagup", "Cup"
        end 
        
        corr = correlation_matrix(ψ, op1, op2; sites=L_end(sys):R_begin(sys))
        return 2 * imag(corr[1, end])
    end 

    workdir = getworkdir(workflag)

    if isnothing(ψ)
        T = []
        currentLR = []

        for file in get_dyna_files()

            ψ = load_ψ(file, workflag)
            curtime = get_time(file)
            append!(T, curtime)

            println("Calculating DPT current, $curtime")

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

function dyna_dptcurrent_mix(; ψ=nothing, sys=DPT_mixed(), workflag = "", kwargs...)

    function work(ψ, ULR, sys)
        # except the DD sites
        corr = correlation_matrix(ψ, "Cdag", "C"; sites=union(L_begin(sys):L_end(sys), R_begin(sys):R_end(sys)))
        
        @assert size(ULR) == size(corr)
        return 2 * imag( sum(ULR .* corr))
    end 


    workdir = getworkdir(workflag)

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

            ψ = load_ψ(file, workflag)
            curtime = get_time(file)
            append!(T, curtime)

            println("Calculating DPT mix current, $curtime")

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

function dyna_SDcurrent(; ψ=nothing, sys :: SD_array=nothing, curtime=nothing, workflag = "", corr_cutoff=Inf, kwargs...)
    function work(ψ, sys)
        
        if  sys.systype == "Fermion" 
            ops = [["Cdag", "C"]]
        else
            ops = [["Cdagup", "Cup"], ["Cdagdn", "Cdn"]]
        end 
        
        corrs = corr_work(ψ, ops, Inf, workflag; corr_cutoff=corr_cutoff)
        current = []

        source = sys.source
        drain = sys.drain
        arr = sys.array

        offset = get_systotal(source) + get_systotal(arr)

        for (ic, corr) in enumerate(corrs)
            if typeof(source) == Reservoir_spatial

                for (couplings..., site) in source.ext_contact
                        

                    coupling = couplings[ic]
                    contact = source.contact
                    @show contact, site
                    currentval = -2  * imag(corr[ contact, site]) * abs(coupling)
                    append!(current, currentval)
                end 

                for (couplings..., site) in drain.ext_contact

                    coupling = couplings[ic]
                    contact = drain.contact + offset
                    @show contact, site
                    currentval = -2 * imag(corr[ site, contact]) * abs(coupling)
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

                #@show sourceinds, draininds, sourceadjinds, drainadjinds

                # the contacts has both the source and drain
                for (couplings..., site) in source.ext_contacts[1]
                    
                    coupling = couplings[ic]
                    currentval = sum(Usource .* corr[ sourceadjinds, site])
                    currentval = -2 * imag(currentval) * abs(coupling)
                    append!(current, currentval)
                end 

                for (couplings..., site) in source.ext_contacts[2]

                    coupling = couplings[ic]
                    currentval = sum(Udrain .* corr[ site, drainadjinds])
                    currentval = -2 * imag(currentval) * abs(coupling)
                    append!(current, currentval)
                end 

            end 
        end 

        # check 
        return current
    end 

    workdir = getworkdir(workflag)

    if isnothing(ψ)
        T = []
        SDcurrents = []

        for file in get_dyna_files()

            ψ = load_ψ(file, workflag)
            curtime = get_time(file)
            append!(T, curtime)

            println("Calculating DPT current, $curtime")

            current = work(ψ, sys)
            @show append!(SDcurrents, [current])

        end 

        writedlm(workdir* "times", T)
        writedlm(workdir* "currentSD", SDcurrents)
    else

        current = work(ψ, sys)
        open( workdir * "currentSD", "a") do io
            writedlm(io, [round.(current, sigdigits=6)])
        end 

    end 


end 



function corr_work(ψ ::MPS,  ops:: Vector{Vector{String}}, curtime, workflag; corr_cutoff=Inf, site_lo = nothing, site_hi = nothing)

    workdir = getworkdir(workflag)
    corrs = []

    for (op1, op2) in ops
        corr = correlation_matrix(ψ, op1, op2)
        #@show corr

        if !isnothing(site_lo) && !isnothing(site_hi)
            corr = corr[site_lo:site_hi, site_lo:site_hi]
        end 

        if curtime < corr_cutoff
            outfile = workdir * "corr" * op1 * op2 * ".h5"

            h5open(outfile, isfile( outfile) ? "r+" : "w") do io

                if !haskey(io, string(curtime))
                    write(io, string(curtime), corr)
                else
                    @warn "duplicate key found, not write"
                end 
                
            end 
        else
            @warn "curtime greater than correlation cutoff time, do not write to file"
        end 
        
        append!(corrs, [corr])
    end 

    return corrs

end 


function dyna_corr(; ψ=nothing, sys=Chain(), curtime=nothing, workflag = "", corr_cutoff=Inf, kwargs...) 

    workdir = getworkdir(workflag)
    ops = ops_determiner(sys)
    site_lo, site_hi = site_determiner(sys)

    if isnothing(ψ)

        T = []
        
        for file in get_dyna_files()

            ψ = load_ψ(file, workflag)
            curtime = get_time(file)

            println("Calculating corr, $curtime")
            _ = corr_work(ψ, ops, curtime, workflag; corr_cutoff=corr_cutoff, site_lo = site_lo, site_hi = site_hi)
            append!(T, curtime)

        end 

        writedlm(workdir* "times", T)

    else

        @show curtime
        _ = corr_work(ψ, ops, curtime, workflag ; corr_cutoff=corr_cutoff, site_lo = site_lo, site_hi = site_hi)

    end 

    return nothing

end 




function dyna_LSRcurrent()

    T = []
    current = []

    workdir = getworkdir(workflag)

    for file in get_dyna_files()

        ψ = load_ψ(file, workflag)
        curtime = get_time(file)
        append!(T, curtime)

        println("Calculating LSR current, $curtime")

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

    workdir = getworkdir(workflag)

    for file in get_dyna_files()

        ψ = load_ψ(file, workflag)
        curtime = get_time(file)
        append!(T, curtime)

        println("Calculating partial contraction, $curtime")

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