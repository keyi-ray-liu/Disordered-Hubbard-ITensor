"""Set parameter dictionary for all future calculations"""
function setpara(;L=22, N=11, int_ee=1.0, int_ne=1.0, t=1.0, ζ=[0.5, 0.5], exch=0.2, 
  decay=0.2, self_nuc=false, disorder=false, sweepdim=500, sweepcnt=50, ex=1, weight=10.0, 
  guess=true, manual=false, itr_dis=[1.0], range=1000, noise=true, method="DMRG", QE=0, xscale=1.0, 
  QN=true, CN=11, QEen=1.0, dp=[], ζ_dp = [], QEoffset = 0.0, output="Default")

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
    "CN" => CN,
    "QEen" => QEen,
    "dp" => dp,
    "ζ_dp" => ζ_dp,
    "QEoffset" => QEoffset,
    "output" => output
  )

  if length(dp) != QE || length(ζ_dp) != QE
    throw(ArgumentError("dp parameter(s) and QE number mismatch"))
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