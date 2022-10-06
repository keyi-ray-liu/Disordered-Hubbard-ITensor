"""Completes one iteration of searches based on the number of excited states. Will return states only when all energies are in ascending order."""
function single_search(para::Dict, disx, disy, lambda)
  # Create N fermion indices

  sweepdim = para["sweepdim"]
  sweepcnt = para["sweepcnt"]
  global_cnt = 60
  ex = para["ex"]
  energies = []
  allvars = []
  vars = []
  states::Vector{MPS} = []
  weight = para["weight"]
  itr_dis = para["itr_dis"]
  allres = []
  noise = para["noise"]
  tol = 1e-8
  
  sites = init_site(para)

  if length(itr_dis) > 1

    # initial state with no disorder
    H = init_ham(para, para["L"], disx.* 0.0, disy.* 0.0, sites)
    psi0 = init_state(para, sites, disx.*0.0, disy.*0.0)

    # iteratively build up the GS guess wavefunction
    for scale in itr_dis
      
      sweeps = Sweeps(sweepcnt)
      setmaxdim!(sweeps, sweepdim)

      if noise
        setnoise!(sweeps, 1E-5)
      end 
      
      setcutoff!(sweeps, 1E-10)

      #_, psi = dmrg(H, psi0, sweeps)
      _, psi = shift_and_invert(H, psi0, sweeps)

      H = init_ham(para, para["L"], disx.* scale, disy.* scale, sites)
      psi0 = psi

    end 
    

  else
    # init with no modification to disorder
    H = init_ham(para, para["L"], disx, disy, sites)
    psi0 = init_state(para, sites, disx, disy)
  end 

  
  # we calculatethe H^dag H as it will be at the LHS of the linear equation
  H2 = contract(H', H; cutoff=1e-12)
  H2 = replaceprime(H2, 2 => 1)


  cur_ex = 1
  cnt = 1

  # changed to while loop for possible retracting operation
  # global cnt to stop excessive looping
  while cur_ex <= ex && cnt <= global_cnt

    cnt += 1
    println("len energy: ", length(energies), "len state: ", length(states), "cur: ", cur_ex)
    println("energies so far", allres)
    println("variances so far", allvars)

    sweeps = Sweeps(sweepcnt)
    setmaxdim!(sweeps, sweepdim)

    if noise
      setnoise!(sweeps, 1E-5)
    end 

    setcutoff!(sweeps, 1E-10)
    #@show sweeps

    # metric is <phi | psi>
    metric = 1
    # Run the DMRG algorithm, returning energy
    while metric >= tol

      print(metric)
      #energy, psi = dmrg(H, psi0, sweeps)
      
      # both H2 and H are passed down. H2 is used fo the LHS, H is used for the RHS
      psi = shift_and_invert(H, H2, psi0, sweeps, lambda; outputlevel=2)
      
      metric = 1 - inner(psi0, psi)
      # new psi is obtained
      psi0 = psi
    
    end 
    
    energy = inner( psi0', H, psi0) / inner( psi0', psi0)
    lambda = energy
    cur_ex += 1



    var = variance(H, psi)

    append!(allres, energy)
    append!(energies, energy)
    append!(states, [psi])
    append!(allvars, var)
    append!(vars, var)
  end 

  return energies, states, allres, allvars, vars
end
