get_time(raw::String) = parse(Float64, SubString(raw, 1 + length(DYNA_STR), length(raw) - length(".h5")))


get_dyna_files() = sort( filter(x->occursin(DYNA_STR,x), readdir(getworkdir())), by=get_time)



""" calculates von neumann and other higher order renyi entropies"""
function entropies(psi::MPS, b::Int, max_order::Int)

    s = siteinds(psi)  
    orthogonalize!(psi, b)
    _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
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

    for s in 2:L - 1
        append!(soi, s)
    end 

    return vectomat([ entropies(ψ, i, max_order) for i in soi]), [ maximum(size(ψ[ j])) for j in soi] 
  
end 



function dyna_EE(; max_order=3)

    T = []
    bonds = []
    SvN = []
    SRenyi = [[] for _ in 2:max_order]
    #sites = []
    workdir = getworkdir()

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

end 

function dyna_occ(; type="Fermion")

    T = []
    
    workdir = getworkdir()

    if type == "Electron"
        occup = []
        occdn = []
    else
        occs =[]
    end 

    for file in get_dyna_files()

        ψ = load_ψ(file)
        t = get_time(file)
        append!(T, t)

        println("Calculating occ, $t")

        if type == "Electron"
            append!(occup, [expect(ψ, "Nup")])
            append!(occdn, [expect(ψ, "Ndn")])
    
        else
            @show occ = expect(ψ, "N")
            append!(occs, [occ])
        end 

    end 

    if type == "Electron"
        writedlm(workdir * "occup", occup)   
        writedlm(workdir * "occdn", occdn)   
  
    else
        writedlm(workdir * "occ", occs)   
    end 

    writedlm(workdir* "times", T)

end 

function dyna_dptcurrent()

    T = []
    currentLR = []

    workdir = getworkdir()

    for file in get_dyna_files()

        ψ = load_ψ(file)
        t = get_time(file)
        append!(T, t)

        println("Calculating DPT current, $t")

        L = div(length(ψ), 2) - 1
        corr = correlation_matrix(ψ, "Cdag", "C"; sites=1:L+1)


        @show append!(currentLR, 2 * imag(corr[L, L + 1]))

        #println(imag(corr))
    end 

    writedlm(workdir* "times", T)
    writedlm(workdir* "currentLR", currentLR)


end 

function dyna_dptcurrent_mix()

    T = []
    currentLR = []

    workdir = getworkdir()

    ks = readdlm( workdir * "ks")
    LR = readdlm( workdir * "LR")

    UL, UR = Uk(1, ks, LR)
    #ULR = UL .* UR'
    ULR = UR .* UL'

    #@show diag(ULR)

    for file in get_dyna_files()

        ψ = load_ψ(file)
        t = get_time(file)
        append!(T, t)

        println("Calculating DPT mix current, $t")

        # except the DD sites
        L = length(ψ) - 2

        corr = correlation_matrix(ψ, "Cdag", "C"; sites=1:L)
        
        @assert size(ULR) == size(corr)
        

        @show current = 2 * imag( sum(ULR .* corr))
        append!(currentLR, current)

    end 

    writedlm(workdir* "times", T)
    writedlm(workdir* "currentLR", currentLR)


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