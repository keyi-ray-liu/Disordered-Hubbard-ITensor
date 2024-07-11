get_time(raw::String) = parse(Float64, SubString(raw, 1 + length(DYNA_STR), length(raw) - length(".h5")))


get_dyna_files() = sort( 
    filter(x-> !occursin("lasttime", x),
    filter(x->occursin(DYNA_STR,x), readdir(getworkdir())))
    , by=get_time)


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



function dyna_occ(; sys=set_Chain(), ψ=nothing, kwargs...)

    workdir = getworkdir()

    if isnothing(ψ)
        if type(sys) == "Electron"
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

            if type(sys)== "Electron"
                append!(occup, [expect(ψ, "Nup")])
                append!(occdn, [expect(ψ, "Ndn")])
        
            else
                @show occ = expect(ψ, "N")
                append!(occs, [occ])
            end 

        end 

        if type(sys) == "Electron"
            writedlm(workdir * "occup", occup)   
            writedlm(workdir * "occdn", occdn)   
    
        else
            writedlm(workdir * "occ", occs)   
        end 

        writedlm(workdir* "times", T)

    else

        if type(sys) == "Fermion"
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

function dyna_dptcurrent(; ψ=nothing, sys=set_DPT(), kwargs...)

    function work(ψ, sys)
        
        if type(sys) == "Fermion" 
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

function dyna_SDcurrent(; ψ=nothing, sys=set_SD(), kwargs...)
    @error "Not completed"
    @assert typeof(sys) == SD_array "sys type wrong!"
    function work(ψ, sys)
        
        if  sys.type == "Fermion" 
            op1, op2 = "Cdag", "C"
        else
            op1, op2 = "Cdagup", "Cup"
        end 
        
        corr = correlation_matrix(ψ, op1, op2; sites=sys.source.contact:sys.s_contact)
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


function dyna_corr(; ψ=nothing, sys=set_Chain(), t=nothing, kwargs...)

    function work(ψ)

        for (op1, op2) in ops
            corr = correlation_matrix(ψ, op1, op2)
            outfile = workdir * "corr" * op1 * op2 * ".h5"
            h5open(outfile, isfile( outfile) ? "r+" : "w") do io
                write(io, string(t), corr)
            end 
        end 

        return nothing

    end 

    if type(sys) == "Fermion"
        ops = [ 
            ["Cdag", "C"],
            ["N", "N"]
        ]

    elseif typeof(sys) == DPT_avg
        ops = [
            ["Cdagup", "Cup"],
            #["Cdagdn", "Cdn"],
            #["Cdagup", "Cdn"],
            #["Cdagdn", "Cup"],
            ["Nup", "Nup"]
        ]
    end 

    workdir = getworkdir()

    if isnothing(ψ)

        T = []
        
        for file in get_dyna_files()

            ψ = load_ψ(file)
            t = get_time(file)
            println("Calculating corr, $t")
            work(ψ)

            append!(T, t)
        end 

        writedlm(workdir* "times", T)

    else
        
        work(ψ)
    end 


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