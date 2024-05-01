gen_state_str(sys::systems) = shuffle([ get_type_dict(type(sys))[i] for i= 1:4 for _ in 1:N(sys)[i] ])

gen_state_str(sys::QE_two) = vcat(QE_str(sys)[1], gen_state_str(sys.chain_only), QE_str(sys)[2])


gen_state_str(sys::DPT) = vcat(shuffle([ get_type_dict(type(sys))[i] for i=1:4 for _ in 1:N(sys)[i]]),  shuffle([ get_type_dict(type(sys))[i] for i=1:4 for _ in 1:N(sys)[i]]), ["Occ", "Emp"])

gen_state_str(sys::DPT_mixed) = gen_state_str(sys.dpt)

gen_state_str(sys::LSR_SIAM) = vcat([ get_type_dict(type(sys))[i] for i=1:4 for _ in 1:N(sys)[i]], ["Emp"], [ get_type_dict(type(sys))[i] for i=1:4 for _ in 1:N(sys)[i]])

gen_state_str(sys::QE_flat_SIAM) = vcat(
    reduce(vcat, [vcat(QE_str(sys)[j], shuffle([ get_type_dict(type(sys))[i] for i= 1:4 for _ in 1:N(sys)[i] ])) for j in 1:legleft(sys)]), 
    ["Emp"], 
    reduce(vcat, [vcat(shuffle([ get_type_dict(type(sys))[i] for i= 1:4 for _ in 1:N(sys)[i] ]), QE_str(sys)[j]) for j in 1 + legleft(sys) : legright(sys) + legleft(sys)])
    )

gen_state_str(sys::QE_parallel) = vcat( gen_state_str(sys.upper),  gen_state_str(sys.lower), ["Emp", "Occ"])

gen_state_str(sys::QE_HOM) = vcat( gen_state_str(sys.upper),  gen_state_str(sys.lower))

gen_state_map(sys::QE_G_SIAM, sites) = Dict( s => gen_state_str(sys.system)[ sitemap(sys)[s] ]  for s in vertices(sites)) 




function gen_state(sys::GQS; QN=true, kwargs...)
    if init(sys) == 1
        return gen_state(sys.chain_only; QN=QN, kwargs...)

    elseif init(sys) == 2
        sites = siteinds(type(sys), get_systotal(sys); conserve_qns=QN)

        @show state_str1 = vcat(["Occ" for _ in 1:N(sys)[2]], ["Emp" for _ in 1:N(sys)[1]])

        @show state_str2 = [ isodd(n) ? "Emp" : "Occ" for n in 1:L(sys)]

        ψ1 = randomMPS(sites, state_str1; linkdims=10)
        ψ2 = randomMPS(sites, state_str2; linkdims=10)

        ψ = sqrt(0.9) * ψ1 + sqrt(0.1) * ψ2
        return ψ

    else
        error("unknonwn init")

    end 
end 

function QE_str(sys::Union{QE_two})

    ex = ["Occ", "Emp"]
    gs = ["Emp", "Occ"]

    if init(sys) == "Left"
        return ex, gs

    elseif init(sys) == "Right"
        return gs, ex

    elseif init(sys) == "Both"
        return ex, ex

    elseif init(sys) == "None"
        return gs, gs

    else
        error("Unrecognized QE init")
    end 

end 

function QE_str(sys::QE_flat_SIAM)

    ex = ["Occ", "Emp"]
    gs = ["Emp", "Occ"]

    if init(sys) == "1"
        return ex, gs, gs, gs

    elseif init(sys) == "2"
        return ex, ex, gs, gs

    elseif init(sys) == "3"
        return ex, ex, ex, gs

    elseif init(sys) == "4"
        return ex, ex, ex, ex

    else
        error("Unrecognized QE init")
    end 

end 


function gen_state(sys::QE_G_SIAM; QN=true, kwargs...)

    g = gen_graph(sys)
    sites = siteinds(type(sys), g; conserve_qns=QN)

    #print(sites)
    #@show length(gen_state_str(sys)), length(get_systotal(sys))
    @show gen_state_map(sys, sites)
    ψ = ttn(sites, g -> gen_state_map(sys, sites)[g])

    return ψ

end 



function gen_state(sys::systems; QN=true, kwargs...)


    sites = siteinds(type(sys), get_systotal(sys); conserve_qns=QN)

    #@show length(gen_state_str(sys)), length(get_systotal(sys))
    @show gen_state_str(sys)
    ψ = randomMPS(sites, gen_state_str(sys) ; linkdims=10)

    return ψ

end 


function gen_hamiltonian(sys::systems)

    res = OpSum()

    res = add_hop!(sys, res)
    res = add_DensityDensity!(sys, res)
    res = add_onsite!(sys, res)
    res = add_qe!(sys, res)
    res = add_HubbardRepulsion!(sys, res)
    res = add_specific_int!(sys, res)
    
    return res
end 


# All the system specific logic is wrapped with multiple-dispatch neighbor functions, we no longer perform these logics within the code





