"""Main function that executes searches on all given cases and writes results to file
Most are self-explanatory, self_nuc determines a constant shift in the final energy, 
guess determines if use a guess initial wf instead of a random one, manual indicates 
if implement the 2D flattening manually
"""
function main(para; kwargs...)

  L = para["L"]
  disorder = para["disorder"]
  disx, disy = setdisorder(disorder, L)
  output = para["output"]
  
  # disable threading in BLAS, so that jobs are multithreaded

  prefix = getworkdir()
  states = get(kwargs, :states, MPS[])
  vars = get(kwargs, :vars, [])

  # in this version, jobs are parallelized and BLAS is not 

  if states == []
    sites = init_site(para; kwargs...)
  else
    sites = siteinds( states[1])
  end

  println("type of sites", typeof(sites))
  println("As of main, type of states", typeof(states))
  
  energy, states, allres, allvars, vars= single_search(para, sites, disx[1, :], disy[1, :], states, vars; kwargs...)

  gaps = energy[2:end] - energy[1:end-1]
  
  writedlm( prefix * output * "ex" , energy)
  writedlm( prefix * output * "gaps", gaps)
  writedlm( prefix * output * "var", vars)
  writedlm( prefix * output * "allvar", allvars)
  writedlm( prefix * output * "allenergy", allres)

  # write wf

  if output == "" 
    wf = h5open( prefix * "wf"  * ".h5", "w")
  else
    wf = h5open( prefix * output * ".h5", "w")
  end 

  for (i, psi) in enumerate(states)
    write(wf, "psi" * string(i), psi)
  end

  close(wf)
  return nothing
end
