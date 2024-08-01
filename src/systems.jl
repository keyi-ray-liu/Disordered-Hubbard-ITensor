const GLOBAL_CNT = 60
const TOL = 1e-8
const DYNA_STR = "tTDVP"
const DPT_INIT_BIAS = [-100.0, 100.0]
const BIAS_LR = 0.5
const LASTSTSTR = "tTDVPlaststate"
const STA_STR = "temp_cur_ex"


abstract type systems end 
abstract type simulations end
abstract type reservoir end

# """ Defines a Hubbard system """
# struct Hubbard <: systems
#     U::Float64
# end 

# U(sys::Hubbard) = sys.U


struct Static <: simulations

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
    
end 

output(simulation::Static) = simulation.output

function set_Static(;
    ex =1,
    prev_state =MPS[], 
    prev_energy =Float64[], 
    prev_var =Float64[], 
    sweepcnt= 60,
    sweepdim =64, 
    noise=true, 
    cutoff=1e-12, 
    krylovdim=8, 
    weight=10,
    output ="wf",
    kwargs...
    )
    
    return Static(
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

struct Dynamic <: simulations

    τ :: Float64
    start :: Float64
    fin :: Float64
    TEcutoff :: Float64
    TEdim :: Int
    nsite :: Int


end 

function set_Dynamic(;
    τ = 0.1,
    start = 0.1,
    fin=20.0,
    TEcutoff=1E-12,
    TEdim=64,
    nsite=2,
    kwargs...
    )

    return Dynamic(
        τ,
        start,
        fin,
        TEcutoff,
        TEdim,
        nsite
    )

end 

""" Defines a Coulomb system """
struct Coulombic <: systems

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

struct Chain_only{T} <: systems where T <: Number

    L :: Int
    N :: Vector{Int}
    type :: String
    t :: Vector{T}
    coulomb :: Coulombic
end 

L(sys::Chain_only) = sys.L
type(sys::Chain_only) = sys.type
N(sys::Chain_only) = sys.N
t(sys::Chain_only) = sys.t
ζ(sys::Chain_only) = sys.coulomb.ζ
get_systotal(sys::Chain_only) = sys.L

function set_Chain(;
    L =2,
    N =1,
    type="Fermion",
    t=-1,
    kwargs...
    )

    t = FermionCondition(type, t)
    coulomb = set_Coulombic(;kwargs...)

    if typeof(N) == Int
        N = [L - N, N, 0, 0]
    end 

    return Chain_only(
        L,
        N,
        type,
        t,
        coulomb
    )

end 

struct biased_chain <: systems

    chain:: Chain_only
    full_size :: Int
    chain_start :: Int

end 

type(sys::biased_chain) = type(sys.chain)
get_systotal(sys::biased_chain) = sys.full_size

set_biased_chain(; chain_start=1, full_size=100, kwargs...) = biased_chain( set_Chain(;kwargs...), full_size, chain_start)

struct GQS <: systems
    chain_only :: Chain_only
    init :: Int
end 

set_GQS(;init=1, kwargs...) = GQS( set_Chain(;kwargs...), init)

init(sys::GQS) = sys.init
N(sys::GQS) = N(sys.chain_only)
L(sys::GQS) = L(sys.chain_only)
type(sys::GQS) = type(sys.chain_only)
get_systotal(sys::GQS) = get_systotal(sys.chain_only)


struct Rectangular{T} <: systems where T <: Number

    Lx :: Int
    Ly :: Int
    N :: Vector{Int}
    type :: String
    t :: Vector{T}
    coulomb :: Coulombic

end 

Lx(sys::Rectangular) = sys.Lx
type(sys::Rectangular) = sys.type
N(sys::Rectangular) = sys.N
t(sys::Rectangular) = sys.t
get_systotal(sys::Rectangular) = sys.Lx * sys.Ly

function set_Rectangular(;
    Lx =3,
    Ly = 3,
    N =1,
    type="Fermion",
    t=-1,
    kwargs...
    )

    t = FermionCondition(type, t)
    coulomb = set_Coulombic(;kwargs...)
    L = Lx * Ly

    if typeof(N) == Int
        N = [L - N, N, 0, 0]
    end 

    return Rectangular(
        Lx,
        Ly,
        N,
        type,
        t,
        coulomb
    )

end 

struct NF_square{T} <: systems where T<: Number

    L :: Int
    N :: Vector{Int}
    U :: Float64
    t :: Vector{T}
    bias :: Float64

end 

get_systotal(sys::NF_square) = sys.L ^2
L(sys::NF_square) = sys.L
t(sys::NF_square) = sys.t
U(sys::NF_square) = sys.U
bias(sys::NF_square) = sys.bias
N(sys::NF_square) = sys.N
type(sys::NF_square) = "Electron"

function set_NF_square(;
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

    return NF_square(
        L,
        [L^2 - Nup - Ndn, Nup, Ndn, 0],
        U,
        t,
        bias
    )


end 


struct DPT <: systems

    U :: Float64
    t_reservoir :: Float64
    t_doubledot :: Float64
    L :: Int
    R :: Int
    type :: String
    contact :: Bool
    contact_t :: Float64
    couple_range :: Int
    bias_doubledot :: Vector{Float64}
    bias_L :: Float64
    bias_R :: Float64
    lattice_info :: Dict

end 


function set_lattice(ddposition, L, R, couple_range)

    lattice_info = Dict{Any, Any}(
        "ddposition" => ddposition
    )

    if ddposition == "R"
        lattice_info["dd_lower"] = L + R + 1
        lattice_info["L_begin"] = 1
        lattice_info["L_end"] = L
        lattice_info["R_begin"] = L + 1
        lattice_info["R_end"] = L + R
        lattice_info["L_contact"] = L - couple_range + 1
        lattice_info["R_contact"] = L + couple_range
        
    elseif ddposition == "M"
        lattice_info["dd_lower"] = L + 1
        lattice_info["L_begin"] = 1
        lattice_info["L_end"] = L
        lattice_info["R_begin"] = L + 3
        lattice_info["R_end"] = L + R + 2
        lattice_info["L_contact"] = L - couple_range + 1
        lattice_info["R_contact"] = L + couple_range + 2

    elseif ddposition == "L"
        lattice_info["dd_lower"] = 1
        lattice_info["L_begin"] = 3
        lattice_info["L_end"] = L + 2
        lattice_info["R_begin"] = L + 3
        lattice_info["R_end"] = L + R + 2
        lattice_info["L_contact"] = L - couple_range + 3
        lattice_info["R_contact"] = L + couple_range + 2

    elseif ddposition == "avg"
        lattice_info["dd_lower"] = L
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


function set_DPT(;
    U = 2.0,
    t_reservoir = 1.0,
    t_doubledot = 1.0,
    L = 6,
    R = 6,
    type = "Fermion",
    couple_range=2,
    contact = true,
    contact_t = 1.0,
    bias_doubledot = [0.0, 0.0], 
    bias_L = 0.0,
    bias_R = 0.0,
    ddposition = "R",
    graph=false,
    kwargs...
    )

    if graph
        ddposition = "M"
    end 

    lattice_info = set_lattice(ddposition, L, R, couple_range)

    sys = DPT(
    U,
    t_reservoir,
    t_doubledot,
    L,
    R,
    type,
    contact,
    contact_t,
    couple_range,
    bias_doubledot,
    bias_L,
    bias_R,
    lattice_info
    )

    if graph
        sys = set_DPT_graph(sys)
    end 

    return sys
end 

type(sys::DPT) = sys.type
U(sys::DPT) = sys.U
L(sys::DPT) = sys.L
R(sys::DPT) = sys.R
couple_range(sys::DPT) = sys.couple_range
bias_doubledot(sys::DPT) = sys.bias_doubledot
bias_L(sys::DPT) = sys.bias_L
bias_R(sys::DPT) = sys.bias_R
t_reservoir(sys::DPT) = sys.t_reservoir
t_doubledot(sys::DPT) = sys.t_doubledot
contact(sys::DPT) = sys.contact
contact_t(sys::DPT) = sys.contact_t
get_systotal(sys::DPT) = 2 + L(sys) + R(sys)
N(sys::DPT) = [L(sys) - div(L(sys), 2), div(L(sys), 2), 0, 0]
dd_lower(sys::DPT) = sys.lattice_info["dd_lower"]
L_begin(sys::DPT) = sys.lattice_info["L_begin"] 
L_end(sys::DPT) = sys.lattice_info["L_end"] 
R_begin(sys::DPT) = sys.lattice_info["R_begin"] 
R_end(sys::DPT) = sys.lattice_info["R_end"] 
L_contact(sys::DPT) = sys.lattice_info["L_contact"]
R_contact(sys::DPT) = sys.lattice_info["R_contact"]
ddposition(sys::DPT) = sys.lattice_info["ddposition"]

struct DPT_mixed <: systems
    dpt :: DPT
    energies :: Vector
    ks :: Vector
    LR :: Vector
    includeU :: Bool
    #U :: Matrix
   # w :: Vector
end 

function set_DPT_mixed(;
    energies = [],
    ks = [],
    LR = [],
    includeU = false,
    ddposition = "R",
    graph = false,
    kwargs...
)

    if graph
        ddposition="M"
        includeU = false
    end 

    dpt = set_DPT(;graph=false, ddposition=ddposition, kwargs...)


    sys= DPT_mixed(
        dpt,
        energies,
        ks,
        LR,
        includeU
       # U,
       # w
    )

    if graph
        sys = set_DPT_graph(sys)
    end 
    
    return sys


end 


get_systotal(sys::DPT_mixed) = get_systotal(sys.dpt)


# this always ties to the L end, regardless of geometry
L_contact(sys::Union{DPT, DPT_mixed}) = L_end(sys) - couple_range(sys) + 1
R_contact(sys::Union{DPT, DPT_mixed}) = R_begin(sys) + couple_range(sys) - 1

includeU(sys::DPT_mixed) = sys.includeU
energies(sys::DPT_mixed) = sys.energies
ks(sys::DPT_mixed) = sys.ks
LR(sys::DPT_mixed) = sys.LR


struct DPT_avg <: systems

    dpt :: Union{DPT, DPT_mixed}

end 


struct DPT_graph <: systems

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

for func ∈ [type, 
    U, L, R, N, couple_range, bias_doubledot, bias_L, bias_R, t_reservoir, t_doubledot, contact, contact_t, dd_lower, L_begin, L_end, R_begin, R_end, L_contact, R_contact, ddposition
    ]
    
    func = Symbol(func)
    @eval $func(sys::DPT_mixed) = $func(sys.dpt)
    @eval $func(sys::DPT_graph) = $func(sys.dpt)
    @eval $func(sys::DPT_avg) = $func(sys.dpt)
end 


for func ∈ [
    includeU, energies, ks, LR
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
    avg :: Bool 
    ;
    ddposition ="R",
    type="Fermion",
    graph = false,
    includeU = false,
    kwargs...
)  
    # for avg, we remove the DD sites and merge it to L end
    if avg
        ddposition = "avg"
        graph = false
        includeU= false
        type = "Electron"
    end 


    if mixed
        sys = set_DPT_mixed(; type=type, ddposition = ddposition, graph=graph, includeU=includeU, kwargs...)

    else
        sys = set_DPT(; type=type, ddposition=ddposition, graph=graph, kwargs...)
    end 

    if avg
        sys = DPT_avg(sys)
    end 

    @show sys

end 



struct reservoir_spatial <: reservoir

    L :: Int
    type :: String
    t :: Vector{Number}
    N :: Vector{Int}
    contact :: Int
    bias :: Float64

end 

struct reservoir_mixed <: reservoir

    energies :: Vector
    type :: String
    ks :: Vector
    LR :: Vector
    bias :: Float64
    N :: Vector{Int}

end 


get_systotal(res::reservoir_spatial) = res.L
get_systotal(res::reservoir_mixed) = length(res.energies)


struct SD_array{T, U, V} <: systems where {T <: reservoir, U <: systems, V <: Any}

    source :: T 
    drain :: T 
    array :: U
    type :: String
    s_contacts :: Array{V}
    d_contacts :: Array{V}

end 

get_systotal(sys::SD_array) = sum( [get_systotal(subsys) for subsys in [sys.source, sys.array, sys.drain]])

for func ∈ [:type
    ]
    @eval $func(sys::SD_array) = sys.$func
end 

function set_reservoir(;
    L = 12,
    t = -1.0,
    type = "Fermion",
    N = 6,
    contact = L,
    bias = 0.0,
    kwargs...)

    t = FermionCondition(type, t)

    if typeof(N) == Int
        N = [L - N, N, 0, 0]
    end 

    reservoir = reservoir_spatial(L, type, t, N, contact, bias)
    return reservoir

end 


function set_SD(
    ;
    L = 12,
    Ns = [0, 0, 0, 1],
    Na = 0,
    Nd = 0,
    contact_scaling = 2.0,
    s_coupling = -0.01,
    d_coupling = -0.01,
    type = "Fermion",
    kwargs...
)

    s_coupling = FermionCondition(type, s_coupling)
    d_coupling = FermionCondition(type, d_coupling)

    s_contacts, d_contacts = set_SD_contacts(s_coupling, d_coupling, contact_scaling, L)


    source = set_reservoir(; L=L, N=Ns, contact = L, type=type,  kwargs...)
    drain = set_reservoir(; L=L, N=Nd, contact = 1, type=type, kwargs...)
    array = set_Rectangular(; N=Na, type=type, kwargs...)

    SD = SD_array(
        source, 
        drain, 
        array, 
        type, 
        s_contacts,
        d_contacts
    )


    return SD

end 




struct LSR_SIAM <: systems

    t_reservoir :: Float64
    t_couple :: Float64
    L :: Int
    R :: Int
    type :: String
    bias_onsite :: Float64
    bias_L :: Float64
    bias_R :: Float64

end 

function set_LSR_SIAM(;
    t_reservoir = 1.0,
    t_couple = 1/sqrt(2),
    L = 32,
    R = 32,
    type = "Fermion",
    bias_onsite = 1.0,
    bias_L = 0.0,
    bias_R = 0.0
    )

    return LSR_SIAM(
    t_reservoir,
    t_couple,
    L,
    R,
    type,
    bias_onsite,
    bias_L,
    bias_R
    )

end 

type(sys::LSR_SIAM) = sys.type
L(sys::LSR_SIAM) = sys.L
R(sys::LSR_SIAM) = sys.R
bias_onsite(sys::LSR_SIAM) = sys.bias_onsite
bias_L(sys::LSR_SIAM) = sys.bias_L
bias_R(sys::LSR_SIAM) = sys.bias_R
t_reservoir(sys::LSR_SIAM) = sys.t_reservoir
t_couple(sys::LSR_SIAM) = sys.t_couple
get_systotal(sys::LSR_SIAM) = 1 + L(sys) + R(sys)

N(sys::LSR_SIAM) = [div(L(sys), 2), div(L(sys), 2), 0, 0]


dis(i::Int, j::Int, sys::systems; range=Inf64) = abs(i- j) <= range ? abs(i - j) : Inf64


function dis(i::Int, j::Int, sys::Rectangular; range=Inf64)

    xi = (i - 1) % Lx(sys) + 1
    yi = (i - 1) ÷ Lx(sys) + 1

    xj = (j - 1) % Lx(sys) + 1
    yj = (j - 1) ÷ Lx(sys) + 1

    vec = norm([xi - xj, yi - yj])
    vec = vec <= range ? vec : Inf64
    
    return vec

end 


CoulombParameters(sys::Union{Chain_only, Rectangular}) = CoulombParameters(sys.coulomb::Coulombic)


CoulombParameters(sys::Coulombic) = sys.λ_ee, sys.λ_ne, sys.exch, sys.scr, sys.range, sys.CN, sys.ζ

SimulationParameters(sys::Static) = sys.ex, sys.prev_state, sys.prev_energy, sys.prev_var, sys.sweepcnt, sys.sweepdim, sys.noise, sys.cutoff,sys.krylovdim, sys.weight

SimulationParameters(sys::Dynamic) = sys.τ, sys.start, sys.fin, sys.TEcutoff, sys.TEdim, sys.nsite

add_specific_int!(sys::systems, res) = res