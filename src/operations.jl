"""Get number of maximum sweeps"""
function getmaxsweep(;L=12, tol=1e-11)
  maxsweepdim = 500

  output = @capture_out begin
    paras = setpara(L=L, N=Int(L//2), sweepdim=maxsweepdim)
    main(paras)
  end

  maxsweep = stringmaxproc(output, tol)
  return maxsweep
end

"""work function that get max sweeps for different lengths of chains"""
function sweepstat()
  res = zeros(0)
  for L in range(12, 100, step=2)
    append!(res, getmaxsweep(L=L))
    println(res)
  end
end

"""work function that runs a series of chains with increasing lengths"""
function plasmonstat(L, lam)

  lam = parse(Float64, lam)
  L = parse(Int, L)
  paras = setpara(L=L, N=div(L, 2), ex=22, int_ee=lam, int_ne=lam, sweepcnt=100, weight=5.0)
  main(paras)


end

"""work function that performs data io and calculates exp vals of observables, in this case, correlation"""
function cal_observe()

  key = "wf"
  workdir = getworkdir()
  obs = "TCD"

  # set range of ex-ex calculations
  cutoff1 = 1
  cutoff2 = 20

  for file in filter(x->occursin(key,x), readdir(workdir))

    wf = h5open( workdir * file, "r")
    exst = sort!(keys(wf), by= x-> parse(Int, x[4:end]))
    numst = length(exst)
    psis = []

    for i = 1:max(cutoff1, cutoff2, numst)
      psi = read(wf, exst[i], MPS)
      append!( psis, [psi])

    end

    for i = 1:min(cutoff1, numst)
      for j = i:min(cutoff2, numst)

        suffix = file[3:4] * "_ex_" * string(i) * "_" *  string(j)
        #cccorr = correlation_matrix(psi,"Cdag","C")

        if obs == "CC"
          obsval = cal_corr(psis[i], psis[j])

        elseif obs == "TCD"
          obsval = cal_tcd(psis[i], psis[j])

        else
          println("unknown observable, abort")
          exit()
        end 

        writedlm( workdir * obs * "_" *  suffix, obsval)

        println("calculation of cc_", suffix, "complete")
      end
    end
  end

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

  L = [2, 30]
  dim = 1000
  cnt = 60
  QN = true

  #QE para
  energy = 0.15
  
  
  QE = 2
  dynamode = "both"
  
  #TE para
  prod = false
  τ = 0.1
  start = 0.1
  fin = 15.0
  TEmethod = "TEBD"
  TEdim = 150
  TEcutoff = 1E-8
  type = "Electron"


  workdir = getworkdir()
  target = "target"

  if !prod && !isfile(workdir * target * ".h5") 
    # if there no target file, we perform a single GS search
    # single search assume no QE GS, headoverride makes sure QE is blocked in Hamiltonian
    paras = setpara(L=L, ex=1, sweepdim=dim, sweepcnt=cnt, QE= QE, QN = QN, output = target, headoverride= (QE > 0) * (QN + 1),
     dynamode=dynamode, type=type)
    main(paras)
  end 
  
  paras = setpara(L=L,  QEen = energy, QE=QE, QEloc=loc, 
  TEcutoff=TEcutoff, TEdim=TEdim, TEmethod=TEmethod, prod=prod, type=type)
  # process wf

  if !prod

    wf = h5open( workdir * target * ".h5", "r")

    if start == 0.1
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

  time_evolve(ψ, sites, paras, start, fin, τ)
  
  
end 


function QE_dynamic_eigen()
  dirs = getworkdir()

end 


function cal_overlap()

  workdir = getworkdir()
  baseline = h5open(workdir * "QE.h5", "r")
  energy = readdlm( workdir * "energy")


  start = 0.1
  fin = 20.0
  step = 0.1

  states = sort(keys(baseline), by= x -> parse(Int, x[4:end]))

  print(states)
  ex = length(states)

  if ex != length(energy)
    error("# of energies do not match # of eigenstates")
  end 

  exwf = Vector{MPS}(undef, ex)

  for (i, key) in enumerate(states)
    exwf[i] = read(baseline, key, MPS)
  end 

  close(baseline)

  L = start:step:fin

  overlap = Array{ComplexF64}(undef, (length(L), ex))
  overlapnorm = zeros((length(L), ex))

  for (i, t) in enumerate(L)

    TE = h5open(workdir * "t" * string(t) *".h5", "r")
    ψ = read( TE, "psi", MPS)
    close(TE)

    for (j, wf) in enumerate(exwf)

      cur = inner(ψ', wf)
      cur *= exp( -im * energy[j] * (t + step))

      overlap[i, j] = cur
      overlapnorm[i, j] = abs(cur)
      println(i, j, " done")

    end 
  end 

  writedlm( workdir * "overlap", overlap)
  writedlm( workdir * "overlapnorm", overlapnorm)

end 

"""calculates the occ for the available t slices"""
function temp_occ(num)
  workdir = getworkdir()
  control = parse(Int, num)

  start = 0.1
  steps = 0.1 
  fin = 10.0

  if control == 1
    method = "TEBD"

  elseif control == 2
    method = "TDVP"

  end 

  res = []
  bonds = []

  T = start:steps:fin

  for t in T
    wf = h5open(workdir * "t$method" * string(t) * ".h5", "r")
    ψ = read( wf, "psi", MPS)
    close(wf)

    bond = checkmaxbond(ψ)
    occ = expect(ψ, "N")
    println(t)
    append!(bonds, bond)
    append!(res, [occ])
  end 

  writedlm(workdir *"time", T)
  writedlm(workdir * "occ", res)
  writedlm(workdir * "bonddim", bonds)
end 


function NF(t, dim, Nup, Ndn)

  t = parse(Float64, t)
  dim = parse(Int, dim)
  Nup = parse(Int, Nup)
  Ndn = parse(Int, Ndn)

  L = [dim, dim]
  Nupdn = 0
  N = [Nup, Ndn, Nupdn]
  ex = 4
  dim = 1000
  cnt = 120
  guess = false
  noise = false

  paras = setpara(L=L, N=N, ex=ex, int_ee=0, int_ne=0, QE=0, t=t,
  guess=guess, sweepdim=dim, sweepcnt=cnt, noise=noise, type="Electron", U=4.0)
  main(paras)

end 