"""Returns a initial state.
Can either generate a random MPS or a guess free fermion initial state"""
function init_state(para, sites, disx, disy)

  L = para["L"]
  N = para["N"]
  guess = para["guess"]
  t = para["t"]
  lambda_ne = para["int_ne"]
  epsilon_ne, epsilon_ee = para["epsilon"]
  decay = para["decay"]
  self = para["self_nuc"]

  # if not guess, using random MPS as init state
  if !guess
    if typeof(L) != Int
      L = L[1] * L[2]
    end 

    # Create an initial random matrix product state
    state = append!([ "Occ" for n=1:N] , ["Emp" for n=1:L -N])
    #@show state
    psi0 = randomMPS(sites,state)


  # use gaussianMPS to calculation an init free fermionic state
  else 
    
    # 1D chain free fermion
    if typeof(L) == Int

      # single particle H
      H = zeros((L, L))

      #hopping
      for i = 1:L-1
        r = dis(i, i + 1, disx, disy)
        hop = hopping(decay, r)
        
        H[i, i + 1] = -t * hop
        H[i + 1, i] = -t * hop
      end 
      
      # nuclear
      for j = 1:L
        if self
          cursum = 0
        else
          cursum = -lambda_ne / epsilon_ne
        end

        for k=1:L
          cursum += lambda_ne / ( dis(j, k, disx, disy) + epsilon_ne)
        end

        H[j, j] = -cursum

      end 

    _, u = eigen(H)

    # Get the Slater determinant
    phi = u[:, 1:N]

    # Create an mps for the free fermion ground state
    psi0 = slater_determinant_to_mps(sites, phi)

    end 


  end 

  return psi0
end
