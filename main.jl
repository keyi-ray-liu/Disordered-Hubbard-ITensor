"""Main function that executes searches on all given cases and writes results to file
Most are self-explanatory, self_nuc determines a constant shift in the final energy, 
guess determines if use a guess initial wf instead of a random one, manual indicates 
if implement the 2D flattening manually
"""
function main(;L=22, N=11, int_ee=1.0, int_ne=1.0, t=1.0, epsilon=[0.5, 0.5], exch=0.2, 
  decay=0.2, self_nuc=false, disorder=false, sweepdim=500, sweepcnt=50, ex=1, weight=10.0, 
  guess=true, manual=false, itr_dis=[1.0], range=1000, noise=true)

  para = setpara(L, N, int_ee, int_ne, t, epsilon, exch, decay, self_nuc, disorder, 
  sweepdim, sweepcnt, ex, weight, guess, manual, itr_dis, range, noise)

  disx, disy = setdisorder(para["disorder"], para["L"])

  # disable threading in BLAS, so that jobs are multithreaded

  prefix = getworkdir()
  cases = size(disx)[1]
  idfile = open( prefix * "disid", "a")

  # in this version, jobs are parallelized and BLAS is not 
  for case=1:cases
    lambda = - 100.0
    energy, states, allres, allvars, vars = single_search(para, disx[case, :], disy[case, :], lambda)

    if typeof(L) != Int
      suffix = string(L[1]) * "x" * string(L[2]) * "_" * string(case)
    else
      suffix = string(L) * "_" * string(case)
    end 
    

    writedlm( prefix * "ex" *  suffix, energy)
    writedlm( prefix * "var" *  suffix, vars)
    writedlm( prefix * "allvar" *  suffix, allvars)
    writedlm( prefix * "allenergy" *  suffix, allres)
    writedlm( idfile, case)

    # write wf
    wf = h5open( prefix * "wf" * suffix * ".h5", "w")

    for (i, psi) in enumerate(states)
      write(wf, "psi" * string(i), psi)
    end

    close(wf)

  end

  close(idfile)
  return nothing
end
