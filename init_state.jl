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

  
  # if not guess, using random MPS as init state
  if !guess
    if typeof(L) != Int
      L = L[1] * L[2]
    end 

    # Create an initial random matrix product state
    QE1 = QE > 0 ? [ "Emp"] : []
    state = append!([ "Occ" for n=1:N] , ["Emp" for n=1:L -N])
    QE2 = QE > 1 ? [ "Emp"] : []

    state = vcat(QE1, state)
    state = vcat(state, QE2)
    #@show state
    ψ0 = randomMPS(sites,state)

    
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

  sites = Vector{Index}(undef, L + QE)

  #left QE
  if QE > 0
    sites[1] = siteind("Fermion"; addtags="n=1", conserve_qns =false)
  end 

  for mid = (QE>0) + 1: (QE>0) + L
    sites[mid] =  siteind("Fermion"; addtags="n=$mid", conserve_qns =QN)
  end 

  if QE > 1
    sites[end] =  siteind("Fermion"; addtags="n=$(L+QE)", conserve_qns =false)
  end 

  #println(length(sites))
  return sites
end 