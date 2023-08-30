

"""work function that runs a series of chains with increasing lengths"""
function plasmonstat(L, lam)

  lam = parse(Float64, lam)
  L = parse(Int, L)
  paras = setpara(L=L, ex=22, int_ee=lam, int_ne=lam, sweepcnt=100, weight=5.0)
  main(paras;)


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

  paras = setpara(L=L, ex=22, int_ee=lam, int_ne=lam, sweepcnt=100, disorder=true, weight=5.0, range=5)
  main(paras;)
end


function GSGap(additional_para)

  paras = setpara(;additional_para..., headoverride=0, QE=0, output="GSGap")
  main(paras;)

end 

"""Statically solve the QE system"""
function QE(QE_para; check_flag = false)

  paras = setpara(;QE_para..., headoverride=0, output="QE")

  if check_flag
    states, energies, vars = load_qe()
    num_state = length(states)
    
    if length(states) < paras["ex"]
      println("number of existing excited states $num_state smaller than required, starting QE")
      main(paras; states=states, energies=energies, vars = vars)

    else
      println("number of existing excited states $num_state already satisfies requirement")
    end 

  else
    main(paras;)
  end 

end 

"""Calculating QE dynamics using various methods"""
function QE_dynamic(simu_para, additional_para)

  #TE para
  product_state = additional_para["product_state"]
  workdir = getworkdir()
  output = "target"

  energy = simu_para[:QEen]
  QN = simu_para[:QN]

  τ = additional_para["τ"]
  start = additional_para["start"]
  fin = additional_para["fin"]
  occ_direct = additional_para["occ_direct"]
  
  if isnothing(energy)
    ex = 2
  else
    ex = 1
  end 

  if !product_state && !isfile(workdir * output * ".h5") 
    # if there no target file, we perform a single GS search
    # single search assume no QE GS, headoverride makes sure QE is blocked in Hamiltonian

    paras = setpara(;simu_para..., ex=ex, output = output, headoverride= QN + 1)
    main(paras;)
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

  time_evolve(ψ, sites, paras, start, fin, occ_direct)
  
  
end 


function NF(t, spdim, dim, Nup, Ndn, geometry)

  t = parse(Float64, t)
  spdim = parse(Int, spdim)
  dim = parse(Int, dim)
  Nup = parse(Int, Nup)
  Ndn = parse(Int, Ndn)

  L = [dim for _ in 1:spdim]
  Nupdn = 0
  N = [Nup, Ndn, Nupdn, 1]
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
  main(paras;)

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

function eigensolver(GS_para, QE_internal_para, QE_para, dyna_para, additional_para)

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
    println("QE file already exists, checking if needs more states")
    QE(QE_para; check_flag = true)
  end 

  #Step 3

  if length(glob("TCD*ex*", workdir)) == 0
    println("TCD file not found, start cal")
    cal_observe(key="QE.h5")
  else
    println("TCD files already exist")
  end 


  #Step 4

  QE_dynamic(dyna_para, additional_para)
  gs_occ()

end 


function time_corr_plot(paras)

  workdir = getworkdir()

  op1 = paras["op1"]
  op2 = paras["op2"]
  tag = paras["tag"]

  files = glob("teigen-wf*.h5", workdir)

  if length(files) == 0
    error(ArgumentError("no time evolved wf found"))
  end 

  get_time(x) = parse(Float64, x[ length(workdir) + length("teigen-wf" ) + 1:end-3])
  sort!(files, by = x -> get_time(x))

  for file in files

    println("calculating file", file)
    wf = h5open(file, "r")
    ψ = read(wf, "psi", MPS)

    NN_corr = abs.(correlation_matrix(ψ, op1, op2))
    writedlm( workdir * tag * "_corr" * string(get_time(file)), NN_corr )
  end 


end 


function SD_dynamics(simu_para, sd_hop, additional_para)


  L = simu_para[:L]
  Ltotal = prod(L)
  if length(L) == 1
    sd_loc = [ [-1.0], [Ltotal]]
    source_site = 1
    drain_site = Ltotal
  else
    y = (minimum(L) - 1) / 2
    sd_loc = [ [y, -1.0 ], [y, maximum(L)]]
    source_site = div(L, 2) + 1
    drain_site = Ltotal - div(L, 2)
  end 

  sd_hop["sd_loc"] = sd_loc
  sd_hop["source_site"] = source_site
  sd_hop["drain_site"] = drain_site

  simu_para[:sd_hop] = sd_hop

  QE_dynamic(simu_para, additional_para)



end 