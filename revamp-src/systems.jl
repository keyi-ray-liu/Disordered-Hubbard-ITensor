
const QESITES = 2
const GLOBAL_CNT = 60
const TOL = 1e-8
const INITIAL_STR = "initialstate"
const DYNA_STR = "tTDVP"
const DPT_INIT_BIAS = [-100.0, 100.0]
const BIAS_LR = 0.25




abstract type systems end 
abstract type simulations end

# """ Defines a Hubbard system """
# struct Hubbard <: systems
#     U::Float64
# end 

# U(sys::Hubbard) = sys.U


struct Static <: simulations

    ex :: Int
    prev_state :: Vector{MPS}
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
    noise=false, 
    cutoff=1e-9, 
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


end 

function set_Dynamic(;
    τ = 0.1,
    start = 0.1,
    fin=20.0,
    TEcutoff=1E-9,
    TEdim=64,
    kwargs...
    )

    return Dynamic(
        τ,
        start,
        fin,
        TEcutoff,
        TEdim
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

struct Chain_only <: systems

    dims :: Int
    N :: Vector{Int}
    type :: String
    t :: Float64
    coulomb :: Coulombic
end 

dims(sys::Chain_only) = sys.dims
type(sys::Chain_only) = sys.type
N(sys::Chain_only) = sys.N
t(sys::Chain_only) = sys.t
ζ(sys::Chain_only) = sys.coulomb.ζ
get_systotal(sys::Chain_only) = sys.dims

function set_Chain(;
    dims =2,
    N =1,
    type="Fermion",
    t=-1,
    coulomb=set_Coulombic(),
    kwargs...
    )

    if typeof(N) == Int
        N = [dims - N, N, 0, 0]
    end 

    return Chain_only(
        dims,
        N,
        type,
        t,
        coulomb
    )

end 

struct Rectangular <: systems

    Lx :: Int
    Ly :: Int
    N :: Vector{Int}
    type :: String
    t :: Float64
    coulomb :: Coulombic

end 

Lx(sys::Rectangular) = sys.dims
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
    coulomb=set_Coulombic(),
    kwargs...
    )

    if typeof(N) == Int
        N = [dims - N, N, 0, 0]
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

struct NF_square <: systems

    L :: Int
    N :: Vector{Int}
    U :: Float64
    t :: Float64
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


    if L^2 - Nup - Ndn != 1
        error("NF electron number incorrect")
    end 

    return NF_square(
        L,
        [1, Nup, Ndn, 0],
        U,
        t,
        bias
    )


end 



"""Regular two QE sys """
struct QE_two <: systems

    QE_distance :: Float64
    offset_scale :: Float64
    QEen :: Float64
    QEmul :: Float64
    dp :: Float64
    init :: String
    product :: Bool
    chain_only :: Chain_only
    
end 

product(sys::QE_two) = sys.product
get_systotal(sys::QE_two) = get_systotal(sys.chain_only) + 2 * QESITES
type(sys::QE_two) = type(sys.chain_only)
t(sys::QE_two) = t(sys.chain_only)
ζ(sys::QE_two) = ζ(sys.chain_only)


function set_QE_two(;
    QE_distance = 2.0,
    offset_scale=0.5,
    QEen=0.0,
    QEmul=1.0,
    dp=1.0,
    init ="Left",
    product=false,
    chain_only=set_Chain(),
    kwargs...
    )

    return QE_two(
        QE_distance,
        offset_scale,
        QEen,
        QEmul,
        dp,
        init,
        product,
        chain_only
    )

end 


"""QE X SIAM system struct, here we assume each leg is attached to each QE, and all legs are connect via a SINGLE 'center' site """
struct QE_flat_SIAM <: systems

    legleft::Int
    legright::Int
    siteseach::Int
    N::Vector{Int}
    type::String
    QE_distance :: Float64
    offset_scale::Float64
    QEen:: Float64
    t :: Float64
    dp :: Float64
    ζ:: Float64
    QEmul :: Float64
    init :: String
    coulomb::Coulombic

end 

legleft(sys::QE_flat_SIAM) = sys.legleft
legright(sys::QE_flat_SIAM) = sys.legright
siteseach(sys::QE_flat_SIAM) = sys.siteseach
N(sys::QE_flat_SIAM) = sys.N
get_systotal(sys::QE_flat_SIAM) = 1 + (legleft(sys) + legright(sys)) * (siteseach(sys) + QESITES)
left(sys::QE_flat_SIAM) = legleft(sys) * (siteseach(sys) + QESITES)
type(sys::QE_flat_SIAM) = sys.type



QE_distance(sys::Union{QE_two, QE_flat_SIAM}) = sys.QE_distance
t(sys::QE_flat_SIAM) = sys.t
QEen(sys::Union{QE_two, QE_flat_SIAM}) = sys.QEen
init(sys::Union{QE_two, QE_flat_SIAM}) = sys.init
QEmul(sys::Union{QE_two, QE_flat_SIAM}) = sys.QEmul
dp(sys::Union{QE_two, QE_flat_SIAM}) = sys.dp
ζ(sys::QE_flat_SIAM) = sys.ζ
offset_scale(sys::Union{QE_two, QE_flat_SIAM}) = sys.offset_scale

function set_QE_SIAM(;
    legleft=2,
    legright=2,
    siteseach=10,
    N=1,
    type="Fermion",
    QE_distance=2,
    offset_scale=0.5,
    QEen=0.0,
    t=-1.0,
    dp=1.0,
    ζ=0.5,
    init="Left",
    QEmul=1.0,
    coulomb=set_Coulombic()
    )

    if typeof(N) == Int
        N = [siteseach - N, N, 0, 0]
    end 

    return QE_flat_SIAM(
        legleft,
        legright,
        siteseach,
        N,
        type,
        QE_distance,
        offset_scale,
        QEen,
        t,
        dp,
        ζ,
        QEmul,
        init,
        coulomb
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

end 

function set_DPT(;
    U = 2.0,
    t_reservoir = 1.0,
    t_doubledot = 1.0,
    L = 32,
    R = 32,
    type = "Fermion",
    couple_range=2,
    contact = true,
    contact_t = 1.0,
    bias_doubledot = [0.0, 0.0], 
    bias_L = 0.0,
    bias_R = 0.0
    )

    return DPT(
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
    bias_R
    )

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

N(sys::DPT) = [div(L(sys), 2), div(L(sys), 2), 0, 0]




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



dis(i, j, sys::systems) = abs(i - j)
dis(i, j, sys::QE_flat_SIAM) = abs(i - j) + QE_distance(sys)
dis(i, j, sys::QE_two) = abs(i - j) + QE_distance(sys)

CoulombParameters(sys::Chain_only) = CoulombParameters(sys.coulomb::Coulombic)

CoulombParameters(sys::QE_flat_SIAM) = CoulombParameters(sys.coulomb::Coulombic)

CoulombParameters(sys::QE_two) = CoulombParameters(sys.coulomb::Coulombic)


CoulombParameters(sys::Coulombic) = sys.λ_ee, sys.λ_ne, sys.exch, sys.scr, sys.range, sys.CN, sys.ζ

SimulationParameters(sys::Static) = sys.ex, sys.prev_state, sys.prev_energy, sys.prev_var, sys.sweepcnt, sys.sweepdim, sys.noise, sys.cutoff,sys.krylovdim, sys.weight

SimulationParameters(sys::Dynamic) = sys.τ, sys.start, sys.fin, sys.TEcutoff, sys.TEdim

add_specific_int!(sys::systems, res) = res