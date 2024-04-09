

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



function NF(paras)

  para = setpara(;paras...)
  main(para;)

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



