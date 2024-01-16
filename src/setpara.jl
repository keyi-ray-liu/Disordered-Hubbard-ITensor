"""Set parameter dictionary for all future calculations, without loading. The system flag now accepts preset parameters defined in here"""
function setpara(;L=22, N="HF", CN="CN", int_ee=2.0, int_ne=2.0, t=-1.0, ζ=[0.5, 0.5], exch=0.2, 
  decay=0.0, self_nuc=false, disorder=false, sweepdim=500, sweepcnt=50, ex=1, weight=10.0, 
  guess=false, manual=false, itr_dis=[1.0], range=10000, noise=true, method="DMRG", QE=0, scales=[1.0], 
  QN=true,  QEen=0.0, dp=1.0, ζ_dp=0.5, QEloc = [], output="", headoverride=0, 
  dynamode="none", cutoff=1E-12, TEcutoff=1E-8, TEdim=300, TEmethod="TDVP", product_state=false, TEBDfactor=2,
  τ=0.1,type="Fermion", U=0.0, snake=false, krylovdim=3, geometry = "linear", spec_hop_struct = Dict{Int64, Float64}(),
  screening_int=0.0, screening_qe=0.0,  s_len = 0, d_len = 0, sd_hop = Dict{Any, Any}(), 
  sd_override=false, range_qe=1000)

  # process L so that it's consistent with the new definition
  if typeof(L) == Int
    L = [L]
  end 

  L = convert(Vector{Int}, L)

  Ltotal = prod(L)
  allnn = get_nn(L, t; snake =snake, geometry=geometry, spec_hop_struct= spec_hop_struct)

  if QE == 0
    dp = []
    ζ_dp = []

  elseif QE == 1
    dp = dp * [1.0]
    ζ_dp = ζ_dp * [1.0]

    if QEloc == []
      if length(L) == 1
        QEloc = [ [-2.0]]
      else
        y = (minimum(L) - 1) / 2
        QEloc = [ [y, -2.0 ]]
      end 
    end 

  elseif QE == 2
    dp = dp * [1.0, 1.0]
    ζ_dp = ζ_dp * [1.0, 1.0]

    if QEloc == []
      if length(L) == 1
        QEloc = [ [-2.0], [prod(L) + 1.0]]
      else
        y = (minimum(L) - 1) / 2
        QEloc = [ [y, -2.0 ], [y, maximum(L) + 1.0]]
      end 
    end 

  end 

  # check half-filling and various N condition
  if N == "HF"

    n = div(prod(L), 2)

    if type == "Electron"
      N = [0, n, n, 0]
    else
      N = [n, n, 0, 0]
    end 

  end



  if typeof(N) == Int
    N = [Ltotal - N, N, 0, 0]
  end 

  if sum(N) > Ltotal
    error("too many sites!")
  end 

  if type == "Fermion" && (N[3] != 0 || N[4] != 0)
    error("Wrong number condition for spinless Fermions!")
  end 

  if length(N) == 3
    append!(N, Ltotal - sum(N))
  end 

  # check charge-neutral
  if CN == "CN"
      CN = N[2] + N[3] + 2 * N[4]
  end 

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
    "snake" => snake,
    "weight" => weight,
    "guess" => guess,
    "manual" => manual,
    "itr_dis" => itr_dis,
    "range" => range,
    "noise" => noise,
    "method" => method,
    "QE" => QE,
    "scales" => scales,
    "QN" => QN,
    "CN" => CN,
    "QEen" => QEen,
    "dp" => dp,
    "ζ_dp" => ζ_dp,
    "QEloc" => QEloc,
    "output" => output,
    "cutoff" => cutoff,
    "headoverride" => headoverride,
    "dynamode" => dynamode,
    "TEcutoff" => TEcutoff,
    "TEdim" => TEdim,
    "TEmethod" => TEmethod,
    "product_state" => product_state,
    "TEBDfactor" => TEBDfactor,
    "τ" => τ,
    "type" => type,
    "U" => U,
    "krylovdim" => krylovdim,
    "screening_int" => screening_int,
    "screening_qe" => screening_qe,
    "s_len" => s_len,
    "d_len" => d_len,
    "sd_hop" => sd_hop,
    "sd_override" => sd_override,
    "range_qe" => range_qe,
    "allnn" => allnn
  )

  if (QE > 0) && s_len + d_len > 0
    error("QE and SD cannot both be non zero")
  end 

  if length(dp) != QE || length(ζ_dp) != QE || length(QEloc) != QE
    error("dp parameter(s) and QE number mismatch")
  end 

  # if dynamode != "none" && headoverride != 2
  #   error("Wrong head override number")
  # end 
  
  workdir = getworkdir()

  open(workdir * output * "parameters.json","w") do f
    JSON3.pretty(f, para)
  end

  #writedlm( workdir * output * "parameters", para)
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