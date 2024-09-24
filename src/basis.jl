gen_mixed(mixed, args...; kwargs...) =  mixed ? gen_mixed(args...;kwargs...) : ([], [], [])

function gen_mixed(L, R, bias_L, bias_R; ordering="SORTED", includeU=true, couple_range=2)
    @info "Set mixed basis, $ordering"
    unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))


    if !includeU
        L -= couple_range
        R -= couple_range
    end 

    # HACK 
    # we test individual version


    L_val = [ (2 *  cos( k * pi / (L + 1) )+ bias_L, k, 1) for k in 1:L] 
    R_val = [ (2 *  cos( k * pi / (R + 1) ) + bias_R, k, -1) for k in 1:R] 


    if ordering == "RNG"

        seed = abs(rand(Int))
        open(getworkdir * "randommixseed", "w") do io
            writedlm(io, seed)
        end 
        result = shuffle( StableRNG(seed), vcat(L_val, R_val))
        
    elseif ordering == "SORTED"
        result = sort( vcat(L_val, R_val), rev=true)

    elseif ordering == "LRSORTED"
        result = vcat(sort(L_val, by= e -> abs(e[1]), rev=true), sort(R_val, by= e -> abs(e[1])))

    elseif ordering == "KMATCHED"
        result = sort( vcat(L_val, R_val), by= e-> e[2])

    elseif ordering == "ABSSORTED"
        result = sort( vcat(L_val, R_val), by = e-> abs(e[1]), rev=true)

    else
        error("Unknown ordering")
    end

    energies, ks, LR= unzip(result)
    writedlm(getworkdir() * "BIASEDenergies", energies)
    writedlm(getworkdir() * "BIASEDenergies" * ordering, energies)

    energies -= LR * bias_L

    writedlm(getworkdir() * "ks", ks)
    writedlm(getworkdir() * "ks" * ordering, ks)
    writedlm(getworkdir() * "LR", LR)
    writedlm(getworkdir() * "LR" * ordering, LR)

    return energies, ks, LR

end 