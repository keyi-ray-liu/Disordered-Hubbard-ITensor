gen_state_str(sys::Systems; kwargs...) = shuffle([ get_type_dict(systype(sys))[i] for i= 1:4 for _ in 1:N(sys)[i] ])

function gen_state_str(sys::QE_two; kwargs...) 


    #state_str = vcat(QE_str(sys)[1], gen_state_str(sys.chain), QE_str(sys)[2])
    state_str = vcat(QE_str(sys)[1], gen_state_str(sys.chain), QE_str(sys)[2])

    return state_str

end 


function gen_state_str(sys::DPT; initdd="LOWER", ifshuffle = false, kwargs...) 

    f(x) = ifshuffle ? shuffle(x) : reverse(x)
    Lres = f([ get_type_dict(systype(sys))[i] for i=1:4 for _ in 1:N(sys)[i]])
    Rres = f([ get_type_dict(systype(sys))[i] for i=1:4 for _ in 1:N(sys)[i]])

    L = []
    M = []
    R = []

    if initdd == "LOWER"
        initarr = ["Occ", "Emp"]
    elseif initdd == "UPPER"
        initarr = ["Emp", "Occ"]
    elseif initdd == "EMPTY"
        initarr = ["Emp", "Emp"]
    else
        error("Unknown init")
    end 

    if ddposition(sys) == "L"
        L = initarr

    elseif ddposition(sys) == "M"
        M = initarr

    else
        R = initarr
    end 

    state_str = vcat(L, Lres, M, Rres, R)
    return state_str
end 


function gen_state_str(sys::DPT_avg; kwargs...)
    Lres = shuffle([ get_type_dict(systype(sys))[i] for i=1:4 for _ in 1:N(sys)[i]])
    Rres = shuffle([ get_type_dict(systype(sys))[i] for i=1:4 for _ in 1:N(sys)[i]])

    if Lres[end] == "Up"
        Lres[end] = "UpDn"

    else
        Lres[end] = "Dn"
    end 


    return vcat(Lres, Rres)

end 

gen_state_str(sys::Reservoir; kwargs...) = reverse([ get_type_dict(sys.systype)[i] for i=1:4 for _ in 1:sys.N[i]])



function gen_state_str(sys::SD_array; fermi=true, random=false, ordering = "SORTED", kwargs...) 

    if fermi && typeof(sys.source) != Reservoir_spatial
        state_str = fermilevel(sys, ordering)

    else
        source = shuffler(gen_state_str(sys.source), random)
        drain = shuffler(gen_state_str(sys.drain), random)

        state_str = vcat( source, gen_state_str(sys.array), drain)
    end 
    
    return state_str
end 



gen_state_str(sys::DPT_mixed; kwargs...) = gen_state_str(sys.dpt; kwargs...)

gen_state_str(sys::LSR_SIAM; kwargs...) = vcat([ get_type_dict(systype(sys))[i] for i=1:4 for _ in 1:N(sys)[i]], ["Emp"], [ get_type_dict(systype(sys))[i] for i=1:4 for _ in 1:N(sys)[i]])

gen_state_str(sys::Biased_chain; kwargs...)  = gen_state_str(sys.chain; kwargs...)


# gen_state_str(sys::QE_flat_SIAM) = vcat(
#     reduce(vcat, [vcat(QE_str(sys)[j], shuffle([ get_type_dict(systype(sys))[i] for i= 1:4 for _ in 1:N(sys)[i] ])) for j in 1:legleft(sys)]), 
#     ["Emp"], 
#     reduce(vcat, [vcat(shuffle([ get_type_dict(systype(sys))[i] for i= 1:4 for _ in 1:N(sys)[i] ]), QE_str(sys)[j]) for j in 1 + legleft(sys) : legright(sys) + legleft(sys)])
#     )

# gen_state_str(sys::QE_parallel) = vcat( gen_state_str(sys.upper),  gen_state_str(sys.lower), ["Emp", "Occ"])

# gen_state_map(sys::QE_G_SIAM, sites) = Dict( s => gen_state_str(sys.system)[ sitemap(sys)[s] ]  for s in vertices(sites)) 

# function QE_str(sys::QE_flat_SIAM)

#     ex = ["Occ", "Emp"]
#     gs = ["Emp", "Occ"]

#     if init(sys) == "1"
#         return ex, gs, gs, gs

#     elseif init(sys) == "2"
#         return ex, ex, gs, gs

#     elseif init(sys) == "3"
#         return ex, ex, ex, gs

#     elseif init(sys) == "4"
#         return ex, ex, ex, ex

#     else
#         error("Unrecognized QE init")
#     end 

# end 

gen_state_str(sys::QE_HOM; kwargs...) = vcat( gen_state_str(sys.upper),  gen_state_str(sys.lower))

gen_state_map(sys::DPT_graph, sites) = Dict( s => gen_state_str(sys.dpt)[ sitemap(sys)[s] ]  for s in vertices(sites)) 


function gen_state(sys::GQS; QN=true, kwargs...)
    if init(sys) == 1
        return gen_state(sys.chain; QN=QN, kwargs...)

    elseif init(sys) == 2
        sites = siteinds(systype(sys), get_systotal(sys); conserve_qns=QN)

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



function gen_state(sys::SD_array; manualmixprod=false, random=false, kwargs...)

    @show manualmixprod
    
    if typeof(sys.source) == Reservoir_spatial || !manualmixprod
        
        state_str =  gen_state_str(sys; random=random, kwargs...)
        @show state_str
        sites = siteinds(systype(sys), get_systotal(sys); conserve_qns=true)

        ψ = randomMPS(sites, state_str
        #; linkdims=10
        )

    else

        
        sites = siteinds(systype(sys), get_systotal(sys); conserve_qns=true)
        source :: Reservoir_momentum = sys.source
        #drain :: Reservoir_momentum = sys.drain
        #array :: Systems = sys.array

        @warn "Only support spatial source unitary for now"

        N = source.N
        Nup = N[2] + N[4]
        Ndn = N[3] + N[4]
        
        #source part LR

        leftind, _ = reservoirmapping(sys)
        @show leftind

        #define the transformation matrices for GaussianMPS
        u = Ujk(source)
        U = zeros( get_systotal(sys), get_systotal(sys))
        U[leftind, leftind] = u

        ϕup = U[:, leftind[1:Nup]]
        ϕdn = U[:, leftind[1:Ndn]]

        ψ = slater_determinant_to_mps(sites, ϕup, ϕdn; maxblocksize=4)
        
        
        @show positiveind(expect(ψ, systype(sys) == "Electron" ? "Ntot" : "N"))


        #drain part LR

        

    end 

    return ψ
end 

# function gen_state(sys::Union{QE_G_SIAM, DPT_graph}; QN=true, kwargs...)
#     @info "gen state for graph"
#     g = gen_graph(sys)
#     sites = siteinds(systype(sys), g; conserve_qns=QN)

#     @show typeof(sites)
#     @show size(sites)
#     #@show length(gen_state_str(sys)), length(get_systotal(sys))

#     #@show gen_state_map(sys, sites)
#     ψ = ttn( g -> gen_state_map(sys, sites)[g], sites)

#     #rng = StableRNG(111)
#     #ψ = normalize(random_ttn(rng, sites; link_space=20))

#     return ψ

# end 




function gen_state(sys::Systems; QN=true, sites=nothing, kwargs...)

    if isnothing(sites)
        sites = siteinds(systype(sys), get_systotal(sys); conserve_qns=QN)

    else
        @info "using predefined sites"
    end 

    #@show length(gen_state_str(sys)), length(get_systotal(sys))
    @show state_str = gen_state_str(sys; kwargs...)
    ψ = randomMPS(sites, state_str
    ; linkdims=10
    )

    #@show expect(ψ, "N")

    return ψ

end 


function gen_hamiltonian(sys::Systems)

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





