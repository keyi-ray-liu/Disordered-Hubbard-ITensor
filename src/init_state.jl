"""Returns a initial state.
Can either generate a random MPS or a guess free fermion initial state"""
function init_state(para, sites, disx, disy)

  L = para["L"]
  N = para["N"]
  guess = para["guess"]
  t = para["t"]
  λ_ne = para["int_ne"]
  ζ_ne, ζ_ee = para["ζ"]
  decay = para["decay"]
  self = para["self_nuc"]
  QE = para["QE"]
  QN = para["QN"]

  
  # if not guess, using random MPS as init state
  if !guess
    if typeof(L) != Int
      L = L[1] * L[2]
    end 

    # Create an initial random matrix product state

    # QN and QE conditions if QN, explicit initiating conditions
    if QN

      AUX1 = QE > 0 ? ["Emp"] : []
      QE1 = QE > 0 ? [ "Occ"] : []
    
      state = append!([ "Occ" for n=1:N] , ["Emp" for n=1:L -N])

      QE2 = QE > 1 ? [ "Occ"] : []
      AUX2 = QE > 1 ? ["Emp"] : []

      state = vcat(QE1, state)
      state = vcat(AUX1, state)
      state = vcat(state, QE2)
      state = vcat(state, AUX2)

      #@show state
      ψ0 = randomMPS(sites,state)

    # random MPS if QN is false, making sure no 'stuck' situation
    else
      ψ0 = randomMPS(sites,L + QE)

    end 

    
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
          cursum = -λ_ne / ζ_ne
        end

        for k=1:L
          cursum += λ_ne / ( dis(j, k, disx, disy) + ζ_ne)
        end

        H[j, j] = -cursum

      end 

      _, u = eigen(H)

    
      # Get the Slater determinant
      ϕ = u[:, 1:N]

      # Create an mps for the free fermion ground state
      ψ0 = slater_determinant_to_mps(sites, ϕ)

    end 


  end 

  return ψ0
end


"""return the sites indices for further use, can add quantum emitters """
function init_site(para::Dict)
  L = para["L"]
  QE = para["QE"]
  QN = para["QN"]

  if typeof(L) != Int
    L = L[1] * L[2]
  end 

  # experimental feature, QN AND QE
  # if we have both, then need AUX sites
  sites = Vector{Index}(undef, L + QE * (QN + 1))

  # determines the number of total sites. IF QE and QN, then we have 2 sites for each QE, else 1
  for s = 1: QE * (QN + 1) + L
    sites[s] =  siteind("Fermion"; addtags="n=$s", conserve_qns =QN)
  end 

  
  #println(length(sites))
  return sites
end 


"""Preparing the state for the TE, if needed"""
function TE_stateprep(ψ, QE)

  if QE == 0
    println("No QE, initial state as read")

  elseif QE == 1
    println(length(ψ))

  elseif QE == 2
    println(length(ψ))

  else
    throw(ArgumentError("invalid QE number"))
  end 

  return ψ

end 