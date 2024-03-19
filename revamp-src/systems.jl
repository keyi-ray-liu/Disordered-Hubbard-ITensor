
const QEsites = 2
abstract type systems end 

""" Defines a Hubbard system """
struct Hubbard <: systems
    U::Float64
end 

U(sys::Hubbard) = sys.U


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

function set_Coulombic(input_json)

    λ_ee = get(input_json, "int_ee", 0.0)
    λ_ne = get(input_json, "int_ne", 0.0)
    exch = get(input_json, "exch", 0.0)
    scr = get(input_json, "scr", 0.0)
    range = get(input_json, "range", 1000)
    CN = get(input_json, "CN", 0.5)
    ζ = get(input_json, "zeta", 0.5)

    return Coulombic(  λ_ee,
    λ_ne,
    exch,
    scr,
    range,
    CN,
    ζ,)

end 

"""Regular two QE sys """
struct QE_two <: systems

    QE_distance :: Float64
    offset_scale :: Float64
    QEen :: Float64
    chain_only :: Chain_only
    
end 


"""QE X SIAM system struct, here we assume each leg is attached to each QE, and all legs are connect via a SINGLE 'center' site """
struct QE_flat_SIAM <: systems

    legleft::Int
    legright::Int
    siteseach::Int
    N::Vector{Int}
    type::String
    offset_scale::Float64
    QEen:: Float64
    coulomb::Coulombic
end 



struct Chain_only <: systems

    dims :: Int
    t :: Float64
    coulomb :: Coulombic
end 

dims(sys::Chain_only) = sys.dims


legleft(sys::QE_flat_SIAM) = sys.legleft
legright(sys::QE_flat_SIAM) = sys.legright
siteseach(sys::QE_flat_SIAM) = sys.siteseach
getsymm(sys::QE_flat_SIAM) = sys.N
type(sys::QE_flat_SIAM) = sys.type

function get_systotal(sys::QE_flat_SIAM)

    return 1 + (legleft(sys) + legright(sys)) * (siteseach(sys) + QEsites)

end 

dis(i, j, sys::QE_flat_SIAM) = abs(i - j)

CoulombParameters(sys::Chain_only) = CoulombParameters(sys.coulomb::Coulombic)

CoulombParameters(sys::QE_flat_SIAM) = CoulombParameters(sys.coulomb::Coulombic)

CoulombParameters(sys::QE_two) = CoulombParameters(sys.coulomb::Coulombic)


CoulombParameters(sys::Coulombic) = sys.λ_ee, sys.λ_ne, sys.exch, sys.scr, sys.range, sys.CN, sys.ζ



add_specifc_int(sys::systems, res) = res