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
function QE_dynamic(simu_para, additional_para; kwargs...)

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

  QEmul = additional_para["QEmul"]
  
  if isnothing(energy)
    ex = 2
  else
    ex = 1
  end 

  if (!product_state || isnothing(energy)) && !isfile(workdir * output * ".h5") 
    # if there no target file, we perform a single GS search
    # single search assume no QE GS, headoverride makes sure QE is blocked in Hamiltonian

    paras = setpara(;simu_para..., ex=ex, output = output, headoverride= QN + 1)
    main(paras; kwargs...)
  end 

  if isnothing(energy)
    plasmon_energy = readdlm(workdir * output * "ex" )
    energy = plasmon_energy[end] - plasmon_energy[end - 1]

    energy = QEmul .* energy


  end 
  
  
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
    paras = setpara(;simu_para...)
    ψ, sites = TE_stateprep(paras)

  end 
  
  print(energy)
  paras = setpara(;simu_para..., τ=τ, QEen=energy, output="TE")
  # further preparation of the initial state, if needed
  #ψ, sites = TE_stateprep(ψ, paras, sites)
  #ψ, sites = TE_stateprep(paras)

  println("length of ψ is:" , length(ψ))
  println("length of sites", length(sites))

  time_evolve(ψ, sites, paras, start, fin, occ_direct; kwargs...)
  
  
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
  
  