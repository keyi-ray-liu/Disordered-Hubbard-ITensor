
shift_basis(left_len, arr_len, vals) = map( x -> x > left_len ? x + arr_len : x , vals)


gen_mixed(mixed :: Bool; kwargs...) =  mixed ? gen_mixed(;kwargs...) : ([], [], [])

function gen_mixed(;L = 4, R =4 , bias_L = 0.0, bias_R=0.0, ω = -1.0,  ordering="SORTED", QPCmixed=true, couple_range=2, workflag = "")
    @info "Set mixed basis, $ordering"
    unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))


    if !QPCmixed
        L -= couple_range
        R -= couple_range
    end 

    # HACK 
    # we test individual version


    L_val = [ (2 *  ω * cos( k * pi / (L + 1) )+ bias_L, k, 1) for k in 1:L] 
    R_val = [ (2 *  ω * cos( k * pi / (R + 1) ) + bias_R, k, -1) for k in 1:R] 


    if ordering == "RNG"

        seed = abs(rand(Int))
        open(getworkdir * "randommixseed", "w") do io
            writedlm(io, seed)
        end 
        result = shuffle( StableRNG(seed), vcat(L_val, R_val))

    elseif ordering == "REVSORTED"
        result = sort( vcat(L_val, R_val), rev=true)

        
    elseif ordering == "SORTED"
        result = sort( vcat(L_val, R_val), rev=false)

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
    writedlm(getworkdir(workflag) * "BIASEDenergies", energies)
    writedlm(getworkdir(workflag) * "BIASEDenergies" * ordering, energies)

    energies -= LR * bias_L

    writedlm(getworkdir(workflag) * "ks", ks)
    writedlm(getworkdir(workflag) * "ks" * ordering, ks)
    writedlm(getworkdir(workflag) * "LR", LR)
    writedlm(getworkdir(workflag) * "LR" * ordering, LR)

    return energies, ks, LR

end 

function Ujk(L::Int, j::Int, k::Int)
    
    v = sqrt(2/ (L + 1)) * sin(j * k * π/ (L + 1))
    return v

end 

function Ujk(sys::Reservoir_momentum, j::Int, k::Int)

    ks = sys.ks
    L = length(ks)
    v = Ujk(L, j, k)

    return v

end 

function Ujk(sys::Reservoir_momentum)

    L = length(sys.ks)
    return vectomat( [[ Ujk(sys, j, k) for j in 1:L] for k in 1:L])

end 





"""generate a product state across LR, representing the 0k fermi level, with all given parameters, at zero bias"""
function fermilevel(sys :: SD_array)

    source = sys.source
    array = sys.array
    drain = sys.drain

    empty = ["Emp" for _ in 1:get_systotal(sys)]

    left_len = get_systotal(source)
    arr_len = get_systotal(array)

    LR = vcat(source.LR, drain.LR)


    # the last occupied sites
    
    sourceinds = findall( x-> x>0, LR)[source.N[1] + 1:end]
    sourceinds = shift_basis(left_len, arr_len, sourceinds)

    draininds = findall( x-> x<0, LR)[drain.N[1] + 1:end]
    draininds = shift_basis(left_len, arr_len, draininds)

    arrayinds = left_len + 1 : left_len + arr_len

    # we assume no mixture of up and dn, for example
    fillstr = systype(sys) == "Fermion" ? "Occ" : "UpDn"

    empty[ draininds ] .= fillstr
    empty[ sourceinds ] .= fillstr

    empty[arrayinds] = gen_state_str(array)

    @show "fermi:", empty
    @show length(empty)

    return empty
    

end 

fermilevel(sys::SD_array, sites)  = random_mps(sites, fermilevel(sys); linkdims=10)
