const GLOBAL_CNT = 400
const TOL = 1e-10
const DYNA_STR = "tTDVP"
const BIASLR = 0.5
const LASTSTSTR = "tTDVPlaststate"
const BACKSTR = "tTDVPbakstate"
const TEMP_tag = "temp_"
const EQINIT_STR = "EqInit"


A = rand()
ITensors.state(::StateName"A", ::SiteType"Fermion") = [sqrt(1 - A), sqrt(A)]


abstract type Subsystem end 
abstract type Systems end 
abstract type SimulationParameters end
abstract type Reservoir end
abstract type ModeDriver end
abstract type TimeControl end
abstract type StateModifier end

# """ Defines a Hubbard system """
# struct Hubbard <: Systems
#     U::Float64
# end 

struct RegularDriver <: ModeDriver
end 

struct BiasGSDriver <: ModeDriver
end 

struct ProdReservoirDriver <: ModeDriver
end 

struct BiasReverseGS <: ModeDriver
end 

struct ProductStateDriver <: ModeDriver
end 


struct Identity <: StateModifier
end 

struct PreLoadGS <: StateModifier
    filestr :: String
end 

struct Supplywf <: StateModifier
    ψ :: MPS
end 


struct LoadSource <: StateModifier
    N :: Int
end 

struct OneStage <: TimeControl
    τ :: Float64
    fin :: Float64
end 

struct TwoStage <: TimeControl
    τ1 :: Float64
    τ2 :: Float64
    fin1 :: Float64
    fin2 :: Float64
end 


# U(sys::Hubbard) = sys.U
systype(sys::Systems) = "Fermion"

struct StaticSimulation <: SimulationParameters

    ex :: Int
    prev_state :: Vector{Any}
    prev_energy :: Vector{Float64}
    prev_var :: Vector{Float64}
    sweepcnt :: Int
    sweepdim :: Int
    noise :: Bool
    cutoff :: Float64
    krylovdim :: Int
    weight :: Float64
    output :: String

    function StaticSimulation(;
        ex =1,
        prev_state =MPS[], 
        prev_energy =Float64[], 
        prev_var =Float64[], 
        sweepcnt= 20,
        sweepdim =64, 
        noise=true, 
        cutoff=1e-12, 
        krylovdim=10, 
        weight=10,
        output ="wf",
        kwargs...
    )
    new(
        ex,
        prev_state,
        prev_energy,
        prev_var,
        sweepcnt,
        sweepdim,
        noise,
        cutoff,
        krylovdim,
        weight,
        output
    )
    end
    
end 

output(simulation::StaticSimulation) = simulation.output


struct DynamicSimulation <: SimulationParameters

    stages :: Array{Array{Float64}}
    # τ :: Float64
    # start :: Float64
    # fin :: Float64
    TEcutoff :: Float64
    TEdim :: Int
    nsite :: Int
    function DynamicSimulation(;  
        τ = 0.1,
        start = 0.0,
        fin=20.0,
        stagetype="uniform",
        TEcutoff=1E-12,
        TEdim=64,
        nsite=2,
        kwargs...)

        if stagetype == "uniform"
            stages = [[τ, start, fin]]

        elseif stagetype == "expadiabatic"
            numstage = get(kwargs, :numstage, 8)
            stepeachstage = get(kwargs, :stepeachstage, 2)
            time = start
            stages = []

            # we first sum up all the residues
            # initstep = τ * (1 - sum( [1/2^num for num in numstage:-1:1]))
            # @show initstep  

            for num in numstage:-1:1
                
                step = τ/2^num

                initstep = num == numstage ? 2 : 1
                duration = step * stepeachstage * initstep
                append!(stages,  [[step, time, time + duration] ])
                time += duration
            end 

            append!(stages, [[τ, time , fin]])
            @show stages
        
        else
            error("Unknown stage control!")
        end 
        new( stages, TEcutoff, TEdim, nsite)
    end 
end 


""" Defines a Coulomb system """
struct Coulombic <: Systems

    λ_ee :: Float64
    λ_ne :: Float64
    exch :: Float64
    scr :: Float64
    range :: Int
    CN :: Float64
    ζ :: Float64
end 

function set_Coulombic(;
    λ_ee = 2.0,
    λ_ne = 2.0,
    exch = 0.2,
    scr =0.0,
    range= 1000,
    CN= 0.5,
    ζ=0.5,
    kwargs...
    )

    return Coulombic(  λ_ee,
    λ_ne,
    exch,
    scr,
    range,
    CN,
    ζ
    )

end 

struct Chain <: Systems

    L :: Int
    N :: Vector{Int}
    systype :: String
    t :: Vector{Number}
    U :: Number
    coulomb :: Coulombic
    function Chain(;
        L =2,
        N =1,
        systype="Fermion",
        t=-1,
        U = 0.0,
        kwargs...
        )
    
        t = FermionCondition(systype, t)
        coulomb = set_Coulombic(;kwargs...)
    
        if typeof(N) == Int
            N = systype == "Electron" ? [L - N, 0, 0, N] : [L - N, N, 0, 0]
        end 
    
        new(
            L,
            N,
            systype,
            t,
            U,
            coulomb
        )

    end 
end 

L(sys::Chain) = sys.L
systype(sys::Chain) = sys.systype
N(sys::Chain) = sys.N
t(sys::Chain) = sys.t
ζ(sys::Chain) = sys.coulomb.ζ
get_systotal(sys::Chain) = sys.L
U(sys::Chain) = sys.U






struct SSH_chain <: Systems

    chain:: Chain
    v :: Float64
    w :: Float64

end 

# for func ∈ [ζ, L, N, systype, get_systotal]
#     @forward SSH_chain.chain func
# end 

ζ(sys::SSH_chain) = ζ(sys.chain)
L(sys::SSH_chain) = L(sys.chain)
N(sys::SSH_chain) = N(sys.chain)
systype(sys::SSH_chain) = systype(sys.chain)
get_systotal(sys::SSH_chain) = get_systotal(sys.chain)
t(sys::SSH_chain, j::Int) = t(sys.chain) .* (isodd(j) ? sys.v : sys.w)

set_SSH_chain(; v=0.0, w=1.0, kwargs...) = SSH_chain( Chain(; kwargs...), v, w)

struct Biased_chain <: Systems

    chain:: Chain
    biaswindow :: Vector{Int}
    bias :: Union{Float64, Int}

end 

systype(sys::Biased_chain) = systype(sys.chain)
get_systotal(sys::Biased_chain) = get_systotal(sys.chain)

set_biased_chain(; biaswindow=[1, 2], bias=-500, kwargs...) = Biased_chain( Chain(;kwargs...), biaswindow, bias)

struct GQS <: Systems
    chain :: Chain
    init :: Int
end 

set_GQS(;init=1, kwargs...) = GQS( Chain(;kwargs...), init)

init(sys::GQS) = sys.init
N(sys::GQS) = N(sys.chain)
L(sys::GQS) = L(sys.chain)
systype(sys::GQS) = systype(sys.chain)
get_systotal(sys::GQS) = get_systotal(sys.chain)


struct Rectangular <: Systems 

    Lx :: Int
    Ly :: Int
    N :: Vector{Int}
    systype :: String
    t :: Vector{Number}
    coulomb :: Coulombic
    U :: Number

    function Rectangular(;
        Lx =3,
        Ly = 3,
        N =1,
        systype="Fermion",
        t=-1,
        U=4.0,
        kwargs...
        )
    
        t = FermionCondition(systype, t)
        coulomb = set_Coulombic(;kwargs...)
        L = Lx * Ly
    
        if typeof(N) == Int
            N = systype == "Electron" ? [L - N, 0, 0, N] : [L - N, N, 0, 0]
        end 

        @show N
    
        new(
            Lx,
            Ly,
            N,
            systype,
            t,
            coulomb,
            U
        )
    
    end 

end 




Lx(sys::Rectangular) = sys.Lx
systype(sys::Rectangular) = sys.systype
N(sys::Rectangular) = sys.N
t(sys::Rectangular) = sys.t
get_systotal(sys::Rectangular) = sys.Lx * sys.Ly
U(sys::Rectangular) = sys.U


struct Ring  <: Systems
    array :: Chain
    function Ring(;
        kwargs...
    )
        array = Chain(;kwargs...)
        new(
            array
        )
    end 

end 


(::Ring)(;kwargs...) = Ring(; kwargs...)

get_systotal(sys::Ring) = get_systotal(sys.array)
N(sys::Ring) = N(sys.array)
U(sys::Ring) = U(sys.array)
systype(sys::Ring) = systype(sys.array)
t(sys::Ring) = t(sys.array)


struct NF_square <: Systems

    L :: Int
    N :: Vector{Int}
    U :: Float64
    t :: Vector{Number}
    bias :: Union{Float64, Int}

    function NF_square(;
        L = 3,
        Nup = 4,
        Ndn = 4,
        U = 4.0,
        t = 0.001,
        bias = 0.0,
        kwargs...
        )
    
        t = FermionCondition("Electron", t)
        if L^2 - Nup - Ndn != 1
            @warn "NF electron (hole) number !=1, is this expected behavior?"
        end 
    
        new(
            L,
            [L^2 - Nup - Ndn, Nup, Ndn, 0],
            U,
            t,
            bias
        )
    
    
    end 

end 

get_systotal(sys::NF_square) = sys.L ^2
L(sys::NF_square) = sys.L
t(sys::NF_square) = sys.t
U(sys::NF_square) = sys.U
bias(sys::NF_square) = sys.bias
N(sys::NF_square) = sys.N
systype(sys::NF_square) = "Electron"



struct DPT <: Systems

    U :: Float64
    t_reservoir :: Float64
    vs :: Float64
    L :: Int
    R :: Int
    systype :: String
    contact :: Bool
    contact_t :: Float64
    couple_range :: Int
    bias_doubledot :: Vector{Float64}
    bias_L :: Float64
    bias_R :: Float64
    lattice_info :: Dict
    μ1 :: Float64
    n1penalty :: Union{Float64, Nothing}

    function DPT(;
    U = 2.0,
    t_reservoir = 1.0,
    vs = 1.0,
    L = 6,
    R = 6,
    systype = "Fermion",
    couple_range=2,
    contact = true,
    contact_t = 1.0,
    bias_doubledot = [0.0, 0.0], 
    bias_L = 0.0,
    bias_R = 0.0,
    ddposition = "R",
    graph=false,
    μ1 = 0.0,
    n1penalty = nothing,
    TLS = false,
    kwargs...
    )

    if graph
        ddposition = "M"
    end 

    lattice_info = set_lattice(ddposition, L, R, couple_range, TLS)
    @show μ1
    new(
    U,
    t_reservoir,
    vs,
    L,
    R,
    systype,
    contact,
    contact_t,
    couple_range,
    bias_doubledot,
    bias_L,
    bias_R,
    lattice_info,
    μ1,
    n1penalty
    )

    # if graph
    #     sys = DPT_graph(sys)
    # end 

end 

end 


function set_lattice(ddposition, L, R, couple_range, TLS)

    lattice_info = Dict{Any, Any}(
        "ddposition" => ddposition
    )

    offset = TLS ? 1 : 0
    if ddposition == "R"
        lattice_info["dd_lower"] = L + R + 1
        lattice_info["dd_upper"] = L + R + 2 - offset
        lattice_info["L_begin"] = 1
        lattice_info["L_end"] = L
        lattice_info["R_begin"] = L + 1
        lattice_info["R_end"] = L + R
        lattice_info["L_contact"] = L - couple_range + 1
        lattice_info["R_contact"] = L + couple_range
        
    elseif ddposition == "M"
        # lattice_info["dd_lower"] = L + 1
        # lattice_info["dd_upper"] = L + 2 - offset
        # lattice_info["L_begin"] = 1
        # lattice_info["L_end"] = L
        # lattice_info["R_begin"] = L + 3 - offset
        # lattice_info["R_end"] = L + R + 2 - offset
        # lattice_info["L_contact"] = L - couple_range + 1
        # lattice_info["R_contact"] = L + couple_range + 2 - offset

        lattice_info["dd_lower"] = L + 1
        lattice_info["dd_upper"] = L + R + 2 - offset
        lattice_info["L_begin"] = 1
        lattice_info["L_end"] = L
        lattice_info["R_begin"] = L + 2 - offset
        lattice_info["R_end"] = L + R + 1 - offset
        lattice_info["L_contact"] = L - couple_range + 1
        lattice_info["R_contact"] = L + couple_range + 1 - offset

    elseif ddposition == "L"
        lattice_info["dd_lower"] = 1
        lattice_info["dd_upper"] = 2 - offset
        lattice_info["L_begin"] = 3 - offset
        lattice_info["L_end"] = L + 2 - offset
        lattice_info["R_begin"] = L + 3 - offset
        lattice_info["R_end"] = L + R + 2 - offset
        lattice_info["L_contact"] = L - couple_range + 3 - offset
        lattice_info["R_contact"] = L + couple_range + 2 - offset

    elseif ddposition == "avg"
        lattice_info["dd_lower"] = L
        lattice_info["dd_upper"] = L + 1
        lattice_info["L_begin"] = 1
        lattice_info["L_end"] = L
        lattice_info["R_begin"] = L + 1
        lattice_info["R_end"] = L + R
        lattice_info["L_contact"] = L - couple_range + 1
        lattice_info["R_contact"] = L + couple_range

    else
        error("Unrecognized dd pos")
    end 

    @show lattice_info
    return lattice_info
end 




systype(sys::DPT) = sys.systype
U(sys::DPT) = sys.U
L(sys::DPT) = sys.L
R(sys::DPT) = sys.R
couple_range(sys::DPT) = sys.couple_range
bias_doubledot(sys::DPT) = sys.bias_doubledot
bias_L(sys::DPT) = sys.bias_L
bias_R(sys::DPT) = sys.bias_R
t_reservoir(sys::DPT) = sys.t_reservoir
vs(sys::DPT) = sys.vs
contact(sys::DPT) = sys.contact
contact_t(sys::DPT) = sys.contact_t
get_systotal(sys::DPT) = 2 + L(sys) + R(sys)
N(sys::DPT) = [L(sys) - div(L(sys), 2), div(L(sys), 2), 0, 0]
dd_lower(sys::DPT) = sys.lattice_info["dd_lower"]
dd_upper(sys::DPT) = sys.lattice_info["dd_upper"]
L_begin(sys::DPT) = sys.lattice_info["L_begin"] 
L_end(sys::DPT) = sys.lattice_info["L_end"] 
R_begin(sys::DPT) = sys.lattice_info["R_begin"] 
R_end(sys::DPT) = sys.lattice_info["R_end"] 
L_contact(sys::DPT) = sys.lattice_info["L_contact"]
R_contact(sys::DPT) = sys.lattice_info["R_contact"]
ddposition(sys::DPT) = sys.lattice_info["ddposition"]

struct DPT_mixed <: Systems
    dpt :: DPT
    energies :: Vector
    ks :: Vector
    LR :: Vector
    QPCmixed :: Bool
    #U :: Matrix
   # w :: Vector
   function DPT_mixed(;
    energies = [],
    ks = [],
    LR = [],
    QPCmixed = false,
    ddposition = "R",
    graph = false,
    kwargs...
)

    if graph
        ddposition="M"
        QPCmixed = false
    end 

    dpt = DPT(;graph=false, ddposition=ddposition, kwargs...)


    new(
        dpt,
        energies,
        ks,
        LR,
        QPCmixed
       # U,
       # w
    )

    # if graph
    #     sys = set_DPT_graph(sys)
    # end 
    
    end 
end 


struct DPT_TLS <: Systems
    dpt :: Union{DPT, DPT_mixed}
end 


struct DPT_TLS2 <: Systems
    dpt :: Union{DPT, DPT_mixed}
end 

get_systotal(sys::Union{DPT_TLS, DPT_TLS2}) = get_systotal(sys.dpt) - 1
get_systotal(sys::DPT_mixed) = get_systotal(sys.dpt)


# this always ties to the L end, regardless of geometry
L_contact(sys::Union{DPT, DPT_mixed}) = L_end(sys) - couple_range(sys) + 1
R_contact(sys::Union{DPT, DPT_mixed}) = R_begin(sys) + couple_range(sys) - 1

QPCmixed(sys::DPT_mixed) = sys.QPCmixed
QPCmixed(sys::DPT_TLS) = QPCmixed(sys.dpt)
QPCmixed(sys::DPT_TLS2) = QPCmixed(sys.dpt)
energies(sys::DPT_mixed) = sys.energies
ks(sys::DPT_mixed) = sys.ks
LR(sys::DPT_mixed) = sys.LR


struct DPT_avg <: Systems

    dpt :: Union{DPT, DPT_mixed}

end 

μ1(sys::DPT_avg) = μ1(sys.dpt)
μ1(sys::DPT_mixed) = μ1(sys.dpt)
μ1(sys::DPT) = sys.μ1


struct DPT_graph <: Systems

    dpt :: Union{DPT, DPT_mixed}
    sitemap :: Dict

end 

function set_DPT_graph(dpt; 
    kwargs...
    )


    ITensors.enable_auto_fermion()
    sitemap = get_sitemap(dpt)
    return DPT_graph(dpt, sitemap)
end 

set_graph(sys::Union{DPT, DPT_mixed}, graph::Bool) = graph ? set_DPT_graph(sys) : sys

for func ∈ [systype, 
    U, L, R, N, couple_range, bias_doubledot, bias_L, bias_R, t_reservoir, vs, contact, contact_t, dd_lower, dd_upper, L_begin, L_end, R_begin, R_end, L_contact, R_contact, ddposition
    ]
    
    func = Symbol(func)
    @eval $func(sys::DPT_mixed) = $func(sys.dpt)
    @eval $func(sys::DPT_graph) = $func(sys.dpt)
    @eval $func(sys::DPT_avg) = $func(sys.dpt)
    @eval $func(sys::DPT_TLS) = $func(sys.dpt)
    @eval $func(sys::DPT_TLS2) = $func(sys.dpt)
end 


for func ∈ [
    QPCmixed, energies, ks, LR
    ]
    
    func = Symbol(func)
    @eval $func(sys::DPT_graph) = $func(sys.dpt::DPT_mixed)
    @eval $func(sys::DPT_avg) = $func(sys.dpt::DPT_mixed)
end 



sitemap(sys::DPT_graph) = sys.sitemap
get_systotal(sys::DPT_graph) = get_systotal(sys.dpt)
get_systotal(sys::DPT_avg) = L(sys) + R(sys)



function DPT_setter(
    mixed :: Bool,
    avg :: Bool,
    TLS :: Bool
    ;
    ddposition ="R",
    systype="Fermion",
    graph = false,
    QPCmixed = false,
    kwargs...
)  
    # for avg, we remove the DD sites and merge it to L end
    if avg
        ddposition = "avg"
        graph = false
        QPCmixed= false
        systype = "Electron"
    end 

    @info ddposition

    if mixed
        sys = DPT_mixed(; systype=systype, ddposition = ddposition, graph=graph, QPCmixed=QPCmixed, TLS = TLS, kwargs...)

    else
        sys = DPT(; systype=systype, ddposition=ddposition, graph=graph, TLS = TLS, kwargs...)
    end 

    if TLS
        sys = DPT_TLS2(sys)
    end 

    if avg
        sys = DPT_avg(sys)
    end 

    #@show sys
    return sys

end 








struct Reservoir_spatial <: Reservoir

    L :: Int
    systype :: String
    t :: Vector{Number}
    N :: Vector{Int}
    contact :: Int
    bias :: Float64
    ext_contact :: Vector

end 

struct Reservoir_momentum <: Reservoir

    energies :: Vector
    systype :: String
    ks :: Vector
    LR :: Vector
    biasS :: Union{Float64, Int}
    biasD :: Union{Float64, Int}
    N :: Vector{Int}
    ext_contacts :: Vector{Vector}

end 


get_systotal(res::Reservoir_spatial) = res.L
get_systotal(res::Reservoir_momentum) = length(res.energies)


struct Plunger <: Subsystem
    onsites :: Array{Number}
    function Plunger( pltype ;G1 = 0.0, G2 = 0.0, Lx = 3, Ly = 3, biasA = 0.0)
        onsites = biasA .* ones(Lx * Ly)

        if pltype == "side"
            onsites[1:Lx] .= G1
            onsites[end - Lx + 1:end] .= G2
        elseif pltype == "diagonal"
            
            if Lx != Ly
                error("Diagonal plunger only supports when Lx = Ly")
            end 
            
            onsites = reshape(onsites, (Lx, Ly))
            for l in 1:Lx - 1

                onsites[ diagind(onsites, l)] .= G1 / (Lx - l)
                onsites[ diagind(onsites, -l)] .= G2 / (Lx - l)
            end 

            onsites = reshape(onsites', Lx * Ly)
        end 
        @show onsites
        new(onsites)

    end 
end 



struct SD_array{T, U} <: Systems where {T <: Reservoir, U <: Systems}

    source :: T 
    drain :: T 
    array :: U
    systype :: String
    plunger :: Plunger

    function SD_array(
        ;
        Ls :: Int= 12,
        Ld :: Int= 12,
        Ns = 1,
        Na = 0,
        Nd  = 0,
        contact_scaling = 2.0,
        s_coupling = -0.01,
        d_coupling = -0.01,
        systype = "Fermion",
        config = "3x3",
        ω = -1.0,  
        G1 = 0.0,
        G2 = 0.0,
        biasA = 0.0,
        pltype = "diagonal",
        kwargs...
        )

        
        s_coupling = FermionCondition(systype, s_coupling)
        d_coupling = FermionCondition(systype, d_coupling)
    
        s_contacts, d_contacts, Lx, Ly = set_SD_parameters(s_coupling, d_coupling, contact_scaling, Ls, config)
    
        @show s_contacts, d_contacts
        @show Ns, Nd
        source, drain = set_reservoir(; Ls=Ls, Ld=Ld, Ns=Ns, Nd=Nd, contacts = [Ls, 1], systype=systype, s_contacts=s_contacts, d_contacts = d_contacts, ω = ω, kwargs...)
    
        if contains(config, "ring")
            #@assert systype == "Electron"
            array = Ring(; Lx=Lx, Ly=Ly, N=Na, systype=systype, kwargs...)
        else
            plunger = Plunger(pltype ; G1 = G1, G2 = G2, Lx = Lx, Ly = Ly, biasA = biasA)
            array = Rectangular(; Lx=Lx, Ly=Ly, N=Na, systype=systype, kwargs...)
        end 
    
        new{Reservoir, Systems}(
            source, 
            drain, 
            array, 
            systype, 
            plunger
        )
    
    end 

end 

get_systotal(sys::SD_array) = sum( [get_systotal(subsys) for subsys in [sys.source, sys.array, sys.drain]])

# for func ∈ [:systype
#     ]
#     @eval $func(sys::SD_array) = sys.$func
# end 

systype(sys::SD_array) = sys.systype

function set_reservoir(;
    Ls :: Int = 12,
    Ld :: Int= Ls,
    Ns = 6,
    Nd = Ns,
    ω = -1.0,
    systype = "Fermion",
    contacts = [Ls, 1],
    reservoir_type = "spatial",
    biasS= 0.0,
    biasD= 0.0,
    energies = [],
    ks = [],
    LR = [],
    s_contacts = [],
    d_contacts = [],
    kwargs...)

    ω = FermionCondition(systype, ω)

    if typeof(Ns) == Int
        Ns = systype == "Electron" ? [Ls - Ns, 0, 0, Ns] : [Ls - Ns, Ns, 0, 0]
    end 

    if typeof(Nd) == Int
        Nd = systype == "Electron" ? [Ld - Nd, 0, 0, Nd] : [Ld - Nd, Nd, 0, 0]
    end 

    #@show Ns, Nd

    if reservoir_type == "spatial"
        source = Reservoir_spatial(Ls, systype, ω, Ns, contacts[1], biasS, s_contacts)
        drain = Reservoir_spatial(Ld, systype, ω, Nd, contacts[2], biasD, d_contacts)

    elseif reservoir_type == "mixed"


        # we assume the reservoir is partitioned according to the Ls, Ld
        source = Reservoir_momentum(energies[1:Ls], systype, ks[1:Ls], LR[1:Ls], biasS, biasD, Ns, [s_contacts, d_contacts])

        drain = Reservoir_momentum(energies[Ls + 1:end], systype, ks[Ls + 1:end], LR[Ls + 1:end], biasS, biasD, Nd, [s_contacts, d_contacts])


    end 
    return source, drain

end 







struct LSR_SIAM <: Systems

    t_reservoir :: Float64
    t_couple :: Float64
    L :: Int
    R :: Int
    systype :: String
    bias_onsite :: Float64
    bias_L :: Float64
    bias_R :: Float64

end 

function set_LSR_SIAM(;
    t_reservoir = 1.0,
    t_couple = 1/sqrt(2),
    L = 32,
    R = 32,
    systype = "Fermion",
    bias_onsite = 1.0,
    bias_L = 0.0,
    bias_R = 0.0
    )

    return LSR_SIAM(
    t_reservoir,
    t_couple,
    L,
    R,
    systype,
    bias_onsite,
    bias_L,
    bias_R
    )

end 

systype(sys::LSR_SIAM) = sys.systype
L(sys::LSR_SIAM) = sys.L
R(sys::LSR_SIAM) = sys.R
bias_onsite(sys::LSR_SIAM) = sys.bias_onsite
bias_L(sys::LSR_SIAM) = sys.bias_L
bias_R(sys::LSR_SIAM) = sys.bias_R
t_reservoir(sys::LSR_SIAM) = sys.t_reservoir
t_couple(sys::LSR_SIAM) = sys.t_couple
get_systotal(sys::LSR_SIAM) = 1 + L(sys) + R(sys)

N(sys::LSR_SIAM) = [div(L(sys), 2), div(L(sys), 2), 0, 0]


dis(i::Int, j::Int, ::Systems; range=Inf64) = abs(i- j) <= range ? abs(i - j) : Inf64


function dis(i::Int, j::Int, sys::Rectangular; range=Inf64)

    xi = (i - 1) % Lx(sys) + 1
    yi = (i - 1) ÷ Lx(sys) + 1

    xj = (j - 1) % Lx(sys) + 1
    yj = (j - 1) ÷ Lx(sys) + 1

    vec = norm([xi - xj, yi - yj])
    vec = vec <= range ? vec : Inf64
    
    return vec

end 


CoulombParameters(sys::Union{Chain, Rectangular}) = CoulombParameters(sys.coulomb::Coulombic)


CoulombParameters(sys::Coulombic) = sys.λ_ee, sys.λ_ne, sys.exch, sys.scr, sys.range, sys.CN, sys.ζ

SimulationParameters(sys::StaticSimulation) = sys.ex, sys.prev_state, sys.prev_energy, sys.prev_var, sys.sweepcnt, sys.sweepdim, sys.noise, sys.cutoff,sys.krylovdim, sys.weight

SimulationParameters(sys::DynamicSimulation) = sys.TEcutoff, sys.TEdim, sys.nsite







# mutable struct SizeObserver <: AbstractObserver
# end

# function ITensors.measure!(o::SizeObserver; bond, sweep, half_sweep, psi, projected_operator, kwargs...)
#   if bond==1 && half_sweep==2
#     psi_size =  Base.format_bytes(Base.summarysize(psi))
#     PH_size =  Base.format_bytes(Base.summarysize(projected_operator))
#     println("After sweep $sweep, |psi| = $psi_size, |PH| = $PH_size")
#   end
# end




function reservoirmapping(sys::SD_array)

    source :: Reservoir_momentum = sys.source
    drain :: Reservoir_momentum = sys.drain
    array :: Systems = sys.array

    offset = get_systotal(source) + get_systotal(array)
    leftind = vcat( positiveind(source.LR), [offset + arr for arr in positiveind(drain.LR)])

    rightind = vcat( negativeind(source.LR), [offset + arr for arr in negativeind(drain.LR)])

    return leftind, rightind


end 


