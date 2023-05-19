

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

"""Statically solve the QE system"""
function QE(num, energy)
  num = parse(Int, num)
  energy = parse(Float64, energy)

  L = 60
  ex = 30
  guess = false
  dim = 300
  cnt = 50
  noise = false


  paras = setpara(L=L, ex=ex, guess=guess, sweepdim=dim, sweepcnt=cnt, noise=noise, QE=num,
  QEen=energy, QEloc=loc, headoverride=0)
  main(paras)

end 

"""Calculating QE dynamics using TEBD"""
function QE_dynamic()

  #L = [2, 30]
  L = 60
  dim = 1000
  cnt = 60
  QN = true

  #QE para
  energy = 0.15672
  dp = 1.0
  
  QE = 2
  dynamode = "both"
  
  #TE para
  product_state = false
  τ = 0.1
  start = 0.1
  fin = 15.0
  TEmethod = "TEBD"
  TEdim = 150
  TEcutoff = 1E-8
  type = "Fermion"
  


  workdir = getworkdir()
  target = "target"

  if !product_state && !isfile(workdir * target * ".h5") 
    # if there no target file, we perform a single GS search
    # single search assume no QE GS, headoverride makes sure QE is blocked in Hamiltonian

    paras = setpara(L=L, ex=1, sweepdim=dim, sweepcnt=cnt, QE= QE, QN = QN, output = target, headoverride= (QE > 0) * (QN + 1),
    dynamode=dynamode, type=type)
    main(paras)
  end 
  
  paras = setpara(L=L,  QEen = energy, QE=QE, TEcutoff=TEcutoff, TEdim=TEdim, dp=dp,
  TEmethod=TEmethod, product_state=product_state, type=type, τ=τ)
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


function NF(t, spdim, dim, Nup, Ndn)

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
  guess=guess, sweepdim=dim, sweepcnt=cnt, noise=noise, type="Electron", U=4.0, snake=snake, krylovdim = krylovdim)
  main(paras)

end 