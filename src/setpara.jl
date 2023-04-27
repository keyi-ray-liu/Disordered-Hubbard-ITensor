"""Set parameter dictionary for all future calculations, without loading. The system flag now accepts preset parameters defined in here"""
function setpara(;L=22, N="HF", CN="CN", int_ee=2.0, int_ne=2.0, t=1.0, ζ=[0.5, 0.5], exch=0.2, 
  decay=0.2, self_nuc=false, disorder=false, sweepdim=500, sweepcnt=50, ex=1, weight=10.0, 
  guess=false, manual=false, itr_dis=[1.0], range=1000, noise=true, method="DMRG", QE=0, scales=[1.0], 
  QN=true,  QEen=1.0, dp=1.0, ζ_dp=0.5, QEloc = [], output="Default", headoverride=0, 
  dynamode="none", TEcutoff=1E-8, TEdim=500, TEmethod="TEBD", product_state=false, TEBDfactor=2,
  τ=0.1,type="Fermion", U=0.0, snake=false, krylovdim=3)

  # process L so that it's consistent with the new definition
  if typeof(L) == Int
    L = [L]
  end 

  if QE == 0
    dp = []
    ζ_dp = []

  elseif QE == 1
    dp = dp * [1.0]
    ζ_dp = ζ_dp * [1.0]

  elseif QE == 2
    dp = dp * [1.0, -1.0]
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

  # check half-filling
  if N == "HF"

    if type == "Electron"
      n = div(prod(L), 2)
      N = [n, n, 0]
    else
      N = div( prod(L), 2)
    end 

  end 

  # check charge-neutral
  if CN == "CN"
    if type == "Electron"
      CN = N[1] + N[2] + 2 * N[3]
    else
      CN = N
    end 
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
    "snake" => snake,
    "krylovdim" => krylovdim
  )

  if (type == "Electron" && typeof(N) == Int) ||  (type == "Fermion" && typeof(N) != Int)
    error("Type must match N: int for fermion, vector for electron")
  end 


  if length(dp) != QE || length(ζ_dp) != QE || length(QEloc) != QE
    error("dp parameter(s) and QE number mismatch")
  end 

  if dynamode != "none" && headoverride != 2
    error("Wrong head override number")
  end 
  
  writedlm("parameters", para)
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