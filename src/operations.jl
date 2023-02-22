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

"""Function to call QE calculations"""
function QE(num, energy)
  num = parse(Int, num)
  energy = parse(Float64, energy)

  L = 60
  N = 30
  CN = 30
  ex = 30
  guess = false
  dim = 300
  cnt = 50
  noise = false
  QN = true
  dps = 1.0
  ζ = 0.5
  offset = 1.0
  
  if num == 1
    dp = dps * [1.0]
    ζ_dp = ζ * [1.0]

  elseif num == 2
    dp = dps * [1.0, -1.0]
    ζ_dp = ζ * [1.0, 1.0]
  end 

  paras = setpara(L=L, N=N, CN=CN, ex=ex, guess=guess, sweepdim=dim, sweepcnt=cnt, noise=noise, QE=num, QN=QN, 
  QEen=energy, dp=dp, ζ_dp=ζ_dp, QEoffset=offset)
  main(paras)

end 

function QE_dynamic()

  L = 12
  N = 6
  CN = 6
  dim = 100
  cnt = 20
  dps = 1.0
  τ = 0.1
  time = 2
  QE = 2
  ζ = 0.5
  energy = 0.5
  offset = 1.0
  QN = true

  if QE == 1
    dp = dps * [1.0]
    ζ_dp = ζ * [1.0]

  elseif QE == 2
    dp = dps * [1.0, -1.0]
    ζ_dp = ζ * [1.0, 1.0]
  end 

  workdir = getworkdir()
  target = "target"

  if !isfile(workdir * target * ".h5")
    # if there no target file, we perform a single GS search
    # single search assume no QE GS
    paras = setpara(L=L, N=N, CN=CN, ex=1, sweepdim=dim, sweepcnt=cnt, output = target, QE=0, QN=QN, headoverride= (QE > 0) * (QN + 1) )
    main(paras)
  end 
  
  paras = setpara(L=L, N=N, CN=CN, sweepdim=dim, sweepcnt=cnt, QEen = energy, QE=QE, QN=QN, dp=dp, ζ_dp=ζ_dp, QEoffset=offset)
  # read the GS
  #wf = h5open( workdir * target * ".h5", "r")
  #ψ = read(wf, "psi1", MPS)

  #sites = siteinds(ψ)
  #ψ = TE_stateprep(ψ, paras, sites)

  # further preparation of the initial state, if needed
  

  ψ, sites = TE_stateprep(paras)
  println("length of ψ is:" , length(ψ))
  println("length of sites", length(sites))

  time_evolve(ψ, sites, paras, time, τ)
  
  
end 