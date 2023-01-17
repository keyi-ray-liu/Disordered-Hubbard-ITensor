"""Set parameter dictionary for all future calculations"""
function setpara(L, N::Int, int_ee::Float64, int_ne::Float64, t::Float64, ζ::Vector{Float64}, exch::Float64, 
  decay::Float64, self_nuc::Bool, disorder::Bool, sweepdim::Int, sweepcnt::Int, ex::Int, weight::Float64, 
  guess::Bool, manual::Bool, itr_dis::Vector{Float64}, range::Int, noise::Bool, method::String, QE::Int, 
  xscale::Float64, QN::Bool, ChargeNeutral::Int, QEen::Float64, dp::Float64)

  # we set the basic parameters for the simulation

  para = Dict(
    "L" => L,
    "N" => N,
    "int_ee" => int_ee,
    "int_ne" => int_ne,
    "t" => t,
    "ζ" => ζ,
    "exch" => exch,
    "decay" => decay,
    "self_nuc" => self_nuc,
    "disorder" => disorder,
    "sweepdim" => sweepdim,
    "sweepcnt" => sweepcnt,
    "ex" => ex,
    "weight" => weight,
    "guess" => guess,
    "manual" => manual,
    "itr_dis" => itr_dis,
    "range" => range,
    "noise" => noise,
    "method" => method,
    "QE" => QE,
    "xscale" => xscale,
    "QN" => QN,
    "CN" => ChargeNeutral,
    "QEen" => QEen,
    "dp" => dp
  )

  if QE > 0 && QN
    throw(ArgumentError("QN and QE condition not compatible"))
  end 
  
  return para
end

"""Set the work directory"""
function getworkdir()

  workdir = pwd() * "/work/"

  if !isdir(workdir)
    mkdir(workdir)
  end 

  return workdir

end 