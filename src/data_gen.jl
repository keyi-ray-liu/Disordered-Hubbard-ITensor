

"""work function that runs a series of chains with increasing lengths"""
function plasmonstat(L, lam)

  lam = parse(Float64, lam)
  L = parse(Int, L)
  paras = setpara(L=L, N=div(L, 2), ex=22, int_ee=lam, int_ne=lam, sweepcnt=100, weight=5.0)
  main(paras)


end



"""work function that either generates disorder or set the parameters for given disorder array"""
function truedisorder(lam, gen)
  lam = parse(Float64, lam)
  gen = parse(Bool, gen)
  workdir = getworkdir()

  if gen
    L = 60
    cases = 24

    disx = (rand(Float64, (cases, L)) .- 0.5) * 2 * 0.3
    disy = (rand(Float64, (cases, L)) .- 0.5) * 2 * 0.2

    writedlm( workdir *  "disx", disx)
    writedlm( workdir * "disy", disy)

  else

    raw = readdlm( workdir * "disx")
    cases, L = size(raw)
  end

  paras = setpara(L=L, N=div(L, 2), ex=22, int_ee=lam, int_ne=lam, sweepcnt=100, disorder=true, weight=5.0, range=5)
  main(paras)
end


function GSGap(additional_para)

  paras = setpara(;additional_para...,
   headoverride=0, ex=2, QE=0, output="GSGap")
  main(paras)

end 

"""Statically solve the QE system"""
function QE(additional_para)

  paras = setpara(;additional_para..., headoverride=0, output="QE")
  main(paras)

end 

"""Calculating QE dynamics using various methods"""
function QE_dynamic(simu_para, time_para)

  #TE para
  product_state = false
  workdir = getworkdir()
  output = "target"

  energy = simu_para[:QEen]
  QN = simu_para[:QN]

  τ = time_para["τ"]
  start = time_para["start"]
  fin = time_para["fin"]
  
  if isnothing(energy)
    ex = 2
  else
    ex = 1
  end 

  if !product_state && !isfile(workdir * output * ".h5") 
    # if there no target file, we perform a single GS search
    # single search assume no QE GS, headoverride makes sure QE is blocked in Hamiltonian

    paras = setpara(;simu_para..., ex=ex, output = output, headoverride= QN + 1)
    main(paras)
  end 

  if isnothing(energy)
    plasmon_energy = readdlm(workdir * output * "ex" )
    energy = plasmon_energy[end] - plasmon_energy[end - 1]
  end 
  
  paras = setpara(;simu_para..., τ=τ, QEen=energy, dynamode="none", output="TE")
  # process wf

  if !product_state

    wf = h5open( workdir * output * ".h5", "r")

    if start == τ
      ψ = read(wf, "psi1", MPS)
    else 
      ψ = read(wf, "psi", MPS)
    end 

    sites = siteinds(ψ)

  else
    ψ, sites = TE_stateprep(paras)

  end 

  # further preparation of the initial state, if needed
  #ψ, sites = TE_stateprep(ψ, paras, sites)
  #ψ, sites = TE_stateprep(paras)

  println("length of ψ is:" , length(ψ))
  println("length of sites", length(sites))

  time_evolve(ψ, sites, paras, start, fin)
  
  
end 


function NF(t, spdim, dim, Nup, Ndn, geometry)

  t = parse(Float64, t)
  spdim = parse(Int, spdim)
  dim = parse(Int, dim)
  Nup = parse(Int, Nup)
  Ndn = parse(Int, Ndn)

  L = [dim for _ in 1:spdim]
  Nupdn = 0
  N = [Nup, Ndn, Nupdn]
  ex = 1
  dim = 1000
  cnt = 120
  guess = false
  noise = false
  snake = true
  krylovdim = 10

  paras = setpara(L=L, N=N, ex=ex, int_ee=0, int_ne=0, QE=0, t=t,
  guess=guess, sweepdim=dim, sweepcnt=cnt, noise=noise, type="Electron", U=4.0, snake=snake, 
  geometry =geometry, krylovdim = krylovdim)
  main(paras)

end 

function eigenplot()

  workdir = getworkdir()

  if !isfile(workdir * "QE.h5")
    raise(ArgumentError("Please provide the eigen basis functions"))
  end 

  staticwffile = h5open( workdir * "QE.h5")
  staticwf=  [read(staticwffile, key, MPS ) for key in sort(keys(staticwffile), by= x-> parse(Int, x[4:end]))]
  println( size(staticwf))

  exps =  [  expect(wf, "N") for wf in staticwf]
  writedlm( workdir * "expects", exps)
  
end 

function eigensolver(GS_para, QE_internal_para, QE_para, dyna_para, time_para)

  workdir = getworkdir()
  
  QE_mul = QE_internal_para["QE_mul"]
  pl_level = QE_internal_para["pl_level"]
  #Step 1

  if !isfile( workdir * "GSGapex")
    println("GS gap file not found, start cal")
    GSGap(GS_para)
  else
    println("GS gap file already exists")
  end 

  plasmon_energy = readdlm(workdir * "GSGapex" )
  QEen = plasmon_energy[pl_level] - plasmon_energy[1]
  QEen *= QE_mul
  #Step 2

  QE_para[:QEen]= QEen 
  dyna_para[:QEen] = QEen

  if !isfile( workdir * "QE.h5")
    println("QE file not found, start cal")
    QE(QE_para)
  else
    println("QE file already exists")
  end 

  #Step 3
  if length(glob("TCD*ex*", workdir)) == 0
    println("TCD file not found, start cal")
    cal_observe(key="QE.h5")
  else
    println("TCD files already exist")
  end 


  #Step 4

  QE_dynamic(dyna_para, time_para)
  gs_occ()

end 
