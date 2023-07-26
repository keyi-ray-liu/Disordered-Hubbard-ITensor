

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


function GSGap()

  L = 120
  N = 6
  guess = false
  dim = 300
  cnt = 50
  noise = false
  krylovdim = 8
  output = "GSGap"

  paras = setpara(L=L, N=N, ex=2, guess=guess, sweepdim=dim, sweepcnt=cnt, noise=noise, QE=0,
   headoverride=0, krylovdim=krylovdim, output=output)
  main(paras)

end 

"""Statically solve the QE system"""
function QE(num, energy)
  num = parse(Int, num)
  energy = parse(Float64, energy)

  L = 120
  N = 6
  ex = 2
  guess = false
  dim = 300
  cnt = 50
  noise = false
  krylovdim = 8

  paras = setpara(L=L, N=N, ex=ex, guess=guess, sweepdim=dim, sweepcnt=cnt, noise=noise, QE=num,
  QEen=energy,  headoverride=0, krylovdim=krylovdim)
  main(paras)

end 

"""Calculating QE dynamics using various methods"""
function QE_dynamic()

  #L = [2, 30]
  L = 400
  dim = 1000
  cnt = 100
  QN = true

  #QE para
  dp = 1.0
  energy = nothing


  QE = 2
  dynamode = "left"
  
  spec_hop_struct = Dict( )
  #TE para
  product_state = false
  τ = 1.0
  start = 1.0
  fin = 100.0
  TEmethod = "TDVP"
  TEdim = 150
  TEcutoff = 1E-8
  type = "Fermion"
  screening = 0.0
  guess = true
  krylovdim = 8
  
  workdir = getworkdir()
  target = "target"
  
  if isnothing(energy)
    ex = 2
  else
    ex = 1
  end 

  if !product_state && !isfile(workdir * target * ".h5") 
    # if there no target file, we perform a single GS search
    # single search assume no QE GS, headoverride makes sure QE is blocked in Hamiltonian

    paras = setpara(L=L, ex=ex, sweepdim=dim, sweepcnt=cnt, QE= QE, QN = QN, output = target, headoverride= (QE > 0) * (QN + 1),
    dynamode=dynamode, type=type, spec_hop_struct=spec_hop_struct, screening=screening, guess=guess, krylovdim=krylovdim)
    main(paras)
  end 

  if isnothing(energy)
    plasmon_energy = readdlm(workdir * "ex1")
    energy = plasmon_energy[end] - plasmon_energy[end - 1]
  end 
  
  paras = setpara(L=L,  QEen = energy, QE=QE, TEcutoff=TEcutoff, TEdim=TEdim, dp=dp, spec_hop_struct=spec_hop_struct, 
  TEmethod=TEmethod, product_state=product_state, type=type, τ=τ, screening=screening)
  # process wf

  if !product_state

    wf = h5open( workdir * target * ".h5", "r")

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