

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



