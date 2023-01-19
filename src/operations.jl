"""Get number of maximum sweeps"""
function getmaxsweep(;L=12, tol=1e-11)
  maxsweepdim = 500

  output = @capture_out begin
    main(L=L, N=Int(L//2), sweepdim=maxsweepdim)
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
  
  main(L=L, N=div(L, 2), ex=22, int_ee=lam, int_ne=lam, sweepcnt=100, weight=5.0)


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

  main(L=L, N=div(L, 2), ex=22, int_ee=lam, int_ne=lam, sweepcnt=100, disorder=true, weight=5.0, range=5)
end

"""Function to call QE calculations"""
function QE(num, energy)
  num = parse(Int, num)
  energy = parse(Float64, energy)

  main(L=12, N=6, CN=6, ex=20, int_ee=2.0, int_ne=2.0, guess=false, method="DMRG", sweepdim=500, sweepcnt=200, noise=false, QE=num, QN=false, QEen=energy, dp=0.0)

end 
