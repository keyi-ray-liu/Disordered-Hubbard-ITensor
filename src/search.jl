"""Completes one iteration of searches based on the number of excited states. Will return states only when all energies are in ascending order."""
function single_search(para::Dict, sites, disx, disy, states, vars; kwargs...)

  prefix = getworkdir()

  energies = get(kwargs, :energies, [])
  λ = get(kwargs, :λ, -20.0)

  sweepdim = para["sweepdim"]
  sweepcnt = para["sweepcnt"]
  global_cnt = 60
  ex = para["ex"]
  allvars = copy(vars)
  weight = para["weight"]
  itr_dis = para["itr_dis"]
  allres = copy(energies)
  noise = para["noise"]
  method = para["method"]
  tol = 1e-8
  krylovdim = para["krylovdim"]
  cutoff = para["cutoff"]
  prefix = getworkdir()



  

  # if we gradually increasing the disorder strength
  if length(itr_dis) > 1

    # initial state with no disorder
    H = init_ham(para, para["L"], disx.* 0.0, disy.* 0.0, sites; kwargs...)
    ϕ = init_state(para, sites, disx.*0.0, disy.*0.0; kwargs...)

    # iteratively build up the GS guess wavefunction
    for scale in itr_dis
      
      sweeps = Sweeps(sweepcnt)
      setmaxdim!(sweeps, sweepdim)

      if noise
        setnoise!(sweeps, 1E-5)
      end 
      
      setcutoff!(sweeps, cutoff)

      _, ψ = dmrg(H, ϕ, sweeps)
      #_, ψ = shift_and_invert(H, ϕ, sweeps)

      H = init_ham(para, para["L"], disx.* scale, disy.* scale, sites; kwargs...)
      ϕ = ψ

    end 


    

  else
    # init with no modification to disorder
    ϕ = init_state(para, sites, disx, disy; kwargs...)
    
    H = init_ham(para, para["L"], disx, disy, sites; kwargs...)

  end 

  
  # we calculatethe H^dag H as it will be at the LHS of the linear equation

  # Oct 12, experiment with official linsolve
  #H2 = contract(H', H; cutoff=1e-12)
  #H2 = replaceprime(H2, 2 => 1)


  cur_ex = 1 + length(states)
  cnt = 1

  # changed to while loop for possible retracting operation
  # global cnt to stop excessive looping
  while cur_ex <= ex && cnt <= global_cnt

    cnt += 1
    println("len energy: ", length(energies), "len state: ", length(states), "cur: ", cur_ex)
    println("energies so far", allres)
    println("variances so far", allvars)

    
    # DMRG block

    if method == "DMRG"

      # DMRG parameters
      sweeps = Sweeps(sweepcnt)
      setmaxdim!(sweeps, sweepdim)

      if noise
        setnoise!(sweeps, 1E-5)
      end 

      setcutoff!(sweeps, cutoff)

      # DMRG method
      # Run the DMRG algorithm, returning energy
      # (dominant eigenvalue) and optimized MPS

      if cur_ex == 1
        energy, ψ = dmrg(H, ϕ, sweeps; eigsolve_krylovdim = krylovdim)
        cur_ex += 1
        

      else
        energy, ψ = dmrg(H, states, ϕ, sweeps; weight=weight, eigsolve_krylovdim = krylovdim) 

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

        # if this is new lowest state, we reset our calculation using ψ as GS
        if cur_ex == 1

          ϕ = ψ
          cur_ex = 2

        end 
        
      end

    
    
      # save temp results
      wf = h5open( prefix * "temp_cur_ex" * string((cur_ex - 1)) * ".h5", "w")
      write(wf, "psi", ψ)
      close(wf)

    # shift and invert block

    elseif method == "SI"
      # metric is <phi | ψ>
      # metric = 1
      itrcnt = 1
      # Run the DMRG algorithm, returning energy

      energy = λ
      # infinity to start as the 'previous' energy
      prev = Inf

      #while metric >= tol

      while itrcnt < 100


        #rewrite shift_and_invert algorithm 
        # now each time, we effectively find (H - λ) |ψ> = |ϕ>
        ψ = shift_and_invert(H, ϕ, λ, sites, para)

        
        #metric = 1 - inner(ϕ, ψ)
        #normalize!(ψ)
        energy = inner( ψ', H, ψ) / inner( ψ', ψ)
        overlap = inner( ψ', ϕ)
        
        println("Power iteration: ", itrcnt, "  Current energy: ", energy, "WF overlap with RHS: ", overlap)
        # new ψ is obtained

        # if the energy converged break the current iteration and start the next
        if abs(prev - energy) <= tol
          break
        end 
        
        ϕ = ψ
        prev = energy
        itrcnt += 1

      end 
      
      cur_ex += 1
      λ = energy



    elseif method == "DMRG-X"

      # tentative application of DMRG X
      energy, ψ = dmrgx(H, ϕ)
      println("energy = ", energy)
      ϕ = ψ

      cur_ex += 1

    else
      println( "Unrecognized method")
      exit()
    end 

    println("As of end of search, type of states", typeof(states))
    var = variance(H, ψ)

    append!(allres, energy)
    append!(energies, energy)
    append!(states, [ψ])
    append!(allvars, var)
    append!(vars, var)
  end 

  return energies, states, allres, allvars, vars
end
