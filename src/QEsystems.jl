
const QESITES = 2
const EMPTY_CENTER = Dict(
    "center_ee" => 0.0,
    "center_ne" => 0.0,
    "center_t" => 0.0,
    "center_range" => 0,
    "center_dis" => 2.0,
    "center_internal_t" => 0.0,
    "center_slope" => 0.5
)

const EMPTY_CONFINE = Dict(
    "confine_start" => 1,
    "confine_range" => Inf,
    "confine_potential" => 0.0
)

const EMPTY_CONFINES = [EMPTY_CONFINE, EMPTY_CONFINE]


"""Regular two QE sys """
struct QE_two <: Systems

    QE_distance :: Float64
    offset_scale :: Float64
    QEen :: Float64
    dp :: Float64
    init :: String
    product :: Bool
    confine_parameter :: Dict
    chain :: Union{Chain, SSH_chain}
    
end 

product(sys::QE_two) = sys.product
get_systotal(sys::QE_two) = get_systotal(sys.chain) + 2 * QESITES
systype(sys::QE_two) = systype(sys.chain)
t(sys::QE_two) = t(sys.chain)
ζ(sys::QE_two) = ζ(sys.chain)
L(sys::QE_two) = L(sys.chain)
N(sys::QE_two) = N(sys.chain)
confine_range(sys::QE_two) = sys.confine_parameter["confine_range"]
confine_potential(sys::QE_two) = sys.confine_parameter["confine_potential"]
confine_start(sys::QE_two) = sys.confine_parameter["confine_start"]

function set_QE_two(;
    QE_distance = 2.0,
    offset_scale=0.5,
    QEen=0.0,
    dp=1.0,
    init ="Left",
    product=false,
    confine_parameters=EMPTY_CONFINES,
    L=10,
    N=5,
    mode ="regular",
    v =0.1,
    w= 1.0,
    kwargs...
    )

    confine_parameter = confine_parameters[1]
    
    if confine_parameter["confine_range"] < L
        N = div(confine_parameter["confine_range"], 2)
        @info "overriding N for confinement, N = $N"
    end 

    if mode == "regular"
        chain = Chain(;L=L, N=N, kwargs...)
    elseif mode == "SSH"
        chain = set_SSH_chain(; L=L, N=N, v=v, w=w, kwargs...)
    end 

    return QE_two(
        QE_distance,
        offset_scale,
        QEen,
        dp,
        init,
        product,
        confine_parameter,
        chain
    )

end 

""" Currently the QE sys 1 and QEsys 2 are put side-by-side , with the center site at the very end"""
struct QE_parallel <: Systems

    upper :: QE_two
    lower :: QE_two
    center_parameter :: Dict

end 

get_upperchain(sys::QE_parallel) = L(sys.upper)
get_lowerchain(sys::QE_parallel) = L(sys.lower)
QEen(sys::QE_parallel) = QEen(sys.upper)
get_systotal(sys::QE_parallel) = get_systotal(sys.upper) + get_systotal(sys.lower) + 2
get_uppertotal(sys::QE_parallel) = get_systotal(sys.upper)
get_lowertotal(sys::QE_parallel) = get_systotal(sys.lower)
systype(sys::QE_parallel) = systype(sys.upper) == systype(sys.lower) ? systype(sys.upper) : error("systype up/lo mismatch!")
offset_scale(sys::QE_parallel) = offset_scale(sys.upper)

function set_QE_parallel(;
    center_parameter = EMPTY_CENTER,
    inits = "1",
    kwargs...
)   

    if inits == "1"
        initstr = ["Left", "None"]
    else
        initstr = ["Left", "Left"]
    end 

    upper = set_QE_two(; init=initstr[1], kwargs...)
    lower = set_QE_two(; init=initstr[2], kwargs...)

    return QE_parallel(
        upper,
        lower,
        center_parameter
    )

end 




"""QE X SIAM system struct, here we assume each leg is attached to each QE, and all legs are connect via a SINGLE 'center' site """
struct QE_flat_SIAM <: Systems

    legleft::Int
    legright::Int
    siteseach::Int
    N::Vector{Int}
    systype::String
    QE_distance :: Float64
    offset_scale::Float64
    QEen:: Float64
    t :: Float64
    dp :: Float64
    ζ:: Float64
    init :: String
    center_parameter :: Dict
    coulomb::Coulombic

end 

legleft(sys::QE_flat_SIAM) = sys.legleft
legright(sys::QE_flat_SIAM) = sys.legright
siteseach(sys::QE_flat_SIAM) = sys.siteseach
N(sys::QE_flat_SIAM) = sys.N
get_systotal(sys::QE_flat_SIAM) = 1 + (legleft(sys) + legright(sys)) * (siteseach(sys) + QESITES)
left(sys::QE_flat_SIAM) = legleft(sys) * (siteseach(sys) + QESITES)
systype(sys::QE_flat_SIAM) = sys.systype


QE_distance(sys::Union{QE_two, QE_flat_SIAM}) = sys.QE_distance
t(sys::QE_flat_SIAM) = sys.t
QEen(sys::Union{QE_two, QE_flat_SIAM}) = sys.QEen
init(sys::Union{QE_two, QE_flat_SIAM}) = sys.init
dp(sys::Union{QE_two, QE_flat_SIAM}) = sys.dp
ζ(sys::QE_flat_SIAM) = sys.ζ
offset_scale(sys::Union{QE_two, QE_flat_SIAM}) = sys.offset_scale

function set_QE_SIAM(;
    legleft=2,
    legright=2,
    siteseach=10,
    N=1,
    systype="Fermion",
    QE_distance=2,
    offset_scale=0.5,
    QEen=0.0,
    t=-1.0,
    dp=1.0,
    ζ=0.5,
    init="1",
    TTN = false,
    center_parameter = EMPTY_CENTER,
    coulomb=set_Coulombic(),
    kwargs...
    )

    if typeof(N) == Int
        N = [siteseach - N, N, 0, 0]
    end 

    sys = QE_flat_SIAM(
        legleft,
        legright,
        siteseach,
        N,
        systype,
        QE_distance,
        offset_scale,
        QEen,
        t,
        dp,
        ζ,
        init,
        center_parameter,
        coulomb
    )

    if !TTN
        return sys
    else

        
        ITensors.enable_auto_fermion()
        sitemap = get_sitemap(sys)
        return QE_G_SIAM(sys, sitemap)
    end

end 

# struct QE_G_SIAM <: Systems
#     system :: QE_flat_SIAM
#     sitemap :: Dict
# end 

# legleft(sys::QE_G_SIAM) = legleft(sys.system)
# legright(sys::QE_G_SIAM) = legright(sys.system)
# siteseach(sys::QE_G_SIAM) = siteseach(sys.system)
# N(sys::QE_G_SIAM) = N(sys.system)
# get_systotal(sys::QE_G_SIAM) = get_systotal(sys.system)
# left(sys::QE_G_SIAM) = left(sys.system)
# systype(sys::QE_G_SIAM) = systype(sys.system)
# t(sys::QE_G_SIAM) = t(sys.system)


# QE_distance(sys::QE_G_SIAM) = QE_distance(sys.system)
# QEen(sys::QE_G_SIAM) = QEen(sys.system)
# init(sys::QE_G_SIAM) = init(sys.system)
# dp(sys::QE_G_SIAM) = dp(sys.system)
# ζ(sys::QE_G_SIAM) = ζ(sys.system)
# offset_scale(sys::QE_G_SIAM) = offset_scale(sys.system)

# sitemap(sys::QE_G_SIAM) = sys.sitemap




struct QE_HOM <: Systems

    upper :: QE_two
    lower :: QE_two
    center_parameter :: Dict
end 

get_upperchain(sys::QE_HOM) = L(sys.upper)
get_lowerchain(sys::QE_HOM) = L(sys.lower)
QEen(sys::QE_HOM) = QEen(sys.upper)
get_systotal(sys::QE_HOM) = get_systotal(sys.upper) + get_systotal(sys.lower) 
get_uppertotal(sys::QE_HOM) = get_systotal(sys.upper)
get_lowertotal(sys::QE_HOM) = get_systotal(sys.lower)
systype(sys::QE_HOM) = systype(sys.upper) == systype(sys.lower) ? systype(sys.upper) : error("systype up/lo mismatch!")
offset_scale(sys::QE_HOM) = offset_scale(sys.upper)

function set_QE_HOM(;
    center_parameter::Dict = EMPTY_CENTER,
    confine_parameters::Vector= EMPTY_CONFINES,
    inits ::String= "1",
    kwargs...
)   

    if inits == "1"
        initstr = ["Left", "None"]
    else
        initstr = ["Left", "Left"]
    end 

    upper = set_QE_two(; init=initstr[1], confine_parameter=confine_parameters[1], kwargs...)
    lower = set_QE_two(; init=initstr[2], confine_parameter=confine_parameters[2], kwargs...)

    return QE_HOM(
        upper,
        lower,
        center_parameter
    )

end 

center_ee(sys::Union{QE_parallel, QE_flat_SIAM, QE_HOM}) = sys.center_parameter["center_ee"]
center_ne(sys::Union{QE_parallel, QE_flat_SIAM, QE_HOM}) = sys.center_parameter["center_ne"]
center_t(sys::Union{QE_parallel, QE_flat_SIAM}) = sys.center_parameter["center_t"]
center_internal_t(sys::QE_parallel) = sys.center_parameter["center_internal_t"]

center_range(sys::Union{QE_parallel, QE_HOM}) = sys.center_parameter["center_range"]
center_dis(sys::Union{QE_parallel, QE_HOM}) = sys.center_parameter["center_dis"]

center_slope(sys::QE_HOM) = sys.center_parameter["center_slope"]

function QE_determiner(key; kwargs...)

    if key == "QE_two"
        sys=  set_QE_two(;  kwargs...)

    elseif key == "QE_SIAM"
        sys = set_QE_SIAM(; kwargs...)

    elseif key == "QE_parallel"
        sys = set_QE_parallel(; kwargs...)

    elseif key == "QE_HOM"
        sys = set_QE_HOM(; kwargs...)

    else
        error("Unrecognized QE key")
    end 

    @show sys
    return sys
end 



function true_center(sys::QE_HOM) 
    uppertotal = get_uppertotal(sys)
    return isodd( uppertotal) ? div(uppertotal, 2) : div(uppertotal, 2) + 1/2
end 

"""this function calculate i, j distance, i in chain 1, j in chain 2. We define an 'X' shape geometry for the y direction"""
function parallel_dis(i, j, sys::QE_HOM) 

    x = abs(i - j)
    y = center_dis(sys) + (abs(i - true_center(sys)) + abs(j - true_center(sys))) * center_slope(sys)

    return sqrt(x^2 + y^2)
end 


qedis(i, j, sys::QE_flat_SIAM) = abs(i -j) + QE_distance(sys)
qedis(i, j, sys::QE_two) = abs(i - j) + QE_distance(sys)


CoulombParameters(sys::QE_flat_SIAM) = CoulombParameters(sys.coulomb::Coulombic)

CoulombParameters(sys::QE_two) = CoulombParameters(sys.coulomb::Coulombic)
