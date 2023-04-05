"""Main function that executes searches on all given cases and writes results to file
Most are self-explanatory, self_nuc determines a constant shift in the final energy, 
guess determines if use a guess initial wf instead of a random one, manual indicates 
if implement the 2D flattening manually
"""
function main(para)

  L = para["L"]
  disorder = para["disorder"]
  disx, disy = setdisorder(disorder, L)
  output = para["output"]
  
  # disable threading in BLAS, so that jobs are multithreaded

  prefix = getworkdir()
  cases = size(disx)[1]
  idfile = open( prefix * "disid", "a")



  # in this version, jobs are parallelized and BLAS is not 
  for case=1:cases
    λ = - 20.0
    sites = init_site(para)

    println("type of sites", typeof(sites))
    energy, states, allres, allvars, vars= single_search(para, sites, disx[case, :], disy[case, :], λ)

    suffix = string(L) * "_" * string(case)
    
    writedlm( prefix * "ex" *  suffix, energy)
    writedlm( prefix * "var" *  suffix, vars)
    writedlm( prefix * "allvar" *  suffix, allvars)
    writedlm( prefix * "allenergy" *  suffix, allres)
    writedlm( idfile, case)

    # write wf

    if output == "Default" || cases > 1
      wf = h5open( prefix * "wf" * suffix * ".h5", "w")
    else
      wf = h5open( prefix * output * ".h5", "w")
    end 

    for (i, psi) in enumerate(states)
      write(wf, "psi" * string(i), psi)
    end

    close(wf)

  end

  close(idfile)
  return nothing
end
