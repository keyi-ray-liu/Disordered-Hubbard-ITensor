get_time(raw::String) = parse(Float64, SubString(raw, 1 + length(DYNA_STR), length(raw) - length(".h5")))


get_dyna_files() = sort( filter(x->occursin(DYNA_STR,x), readdir(getworkdir())), by=get_time)


function entropy_von_neumann(psi::MPS, b::Int)

    s = siteinds(psi)  
    orthogonalize!(psi, b)
    _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
    SvN = 0.0
    for n in 1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
    return SvN
end

function scan_ee(ψ::MPS)

    L = length(ψ)
    soi = []

    for s in 2:L - 1
        append!(soi, s)
    end 
  
    return [ entropy_von_neumann(ψ, i ) for i in soi], [ maximum(size(ψ[ j])) for j in soi], soi
  
end 



function dyna_EE()

    T = []
    bonds = []
    ees = []
    sites = []
    workdir = getworkdir()

    for file in get_dyna_files()

        ψ = load_ψ(file)
        t = get_time(file)

        println("Calculating EE, $t")
        @show ee, bond, site = scan_ee(ψ)

        
        append!(T, t)
        append!(ees, [ee])
        append!(bonds, [bond])
        append!(sites, [site])

    end 

    writedlm(workdir * "times", T)
    writedlm(workdir * "EE", ees)
    writedlm(workdir * "bonds", bonds)
    writedlm(workdir * "sites", sites)

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
        corr = correlation_matrix(ψ, "Cdag", "C")

        @show append!(currentLR, 2 * imag(corr[L, L + 1]))

        #println(imag(corr))
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