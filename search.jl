"""Completes one iteration of searches based on the number of excited states. Will return states only when all energies are in ascending order."""
function single_search(para::Dict, disx, disy)
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
  tol = 1e-6
  
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

      _, psi = dmrg(H, psi0, sweeps)

      H = init_ham(para, para["L"], disx.* scale, disy.* scale, sites)
      psi0 = psi

    end 
    

  else
    # init with no modification to disorder
    H = init_ham(para, para["L"], disx, disy, sites)
    psi0 = init_state(para, sites, disx, disy)
  end 


  cur_ex = 1
  cnt = 1

  # changed to while loop for possible retracting operation
  # global cnt to stop excessive looping
  while cur_ex <= ex && cnt <= global_cnt

    cnt += 1
    println("len energy: ", length(energies), "len state: ", length(states), "cur: ", cur_ex)
    println("energies so far", allres)
    println("variances so far", allvars)
    # Input operator terms which define
    # a Hamiltonian matrix, and convert
    # these terms to an MPO tensor network
    # (1D Hubbard Chain)

    # returns an MPS object
    

    # returns an MPO object


    #@show psi0
    # Plan to do 20 passes or 'sweeps' of DMRG,
    # setting maximum MPS internal dimensions
    # for each sweep and maximum truncation cutoff
    # used when adapting internal dimensions:
    sweeps = Sweeps(sweepcnt)
    setmaxdim!(sweeps, sweepdim)

    if noise
      setnoise!(sweeps, 1E-5)
    end 

    setcutoff!(sweeps, 1E-10)
    #@show sweeps

    # Run the DMRG algorithm, returning energy
    # (dominant eigenvalue) and optimized MPS

    if cur_ex == 1
      energy, psi = dmrg(H, psi0, sweeps)
      cur_ex += 1

    else
      energy, psi = dmrg(H, states, psi0, sweeps; weight) 

      # check if cur energy is lower than previously achieved energy, if so, return to the point with lower energy (if not, start with current state as GS)
      if abs(energy - energies[end]) > tol && energy < energies[end]

        cur = 1

        while cur <= length(energies)

          if abs(energy - energies[cur]) > tol && energies[cur] < energy
            cur += 1
          else 
            break
          end 

        end 

        # reset the current array of states and energies, reset ex count
        energies = energies[begin:cur-1]
        states = states[begin:cur-1]
        vars = vars[begin:cur-1]
        cur_ex = cur
      
      # else continue evaluation
      else
        cur_ex += 1

      end 

    end

    var = variance(H, psi)

    append!(allres, energy)
    append!(energies, energy)
    append!(states, [psi])
    append!(allvars, var)
    append!(vars, var)
  end 

  return energies, states, allres, allvars, vars
end
