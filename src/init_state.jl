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
  headoverride = para["headoverride"]

  
  # if not guess, using random MPS as init state
  if !guess
    if typeof(L) != Int
      L = L[1] * L[2]
    end 

    # Create an initial random matrix product state

    # QN and QE conditions if QN, explicit initiating conditions
    if QN

      state = append!([ "Occ" for n=1:N] , ["Emp" for n=1:L -N]) 

      if headoverride == 0

        AUX1 = QE > 0 ? ["Emp"] : []
        QE1 = QE > 0 ? [ "Occ"] : []
        QE2 = QE > 1 ? [ "Occ"] : []
        AUX2 = QE > 1 ? ["Emp"] : []

        state = vcat(QE1, state)
        state = vcat(AUX1, state)
        state = vcat(state, QE2)
        state = vcat(state, AUX2)

      else

        front = back = ["Emp" for n=1: headoverride]
        state = vcat(front, state)
        state = vcat(state, back)

      end 

      #@show state
      ψ0 = randomMPS(sites,state)

    # random MPS if QN is false, making sure no 'stuck' situation
    else

      ψ0 = randomMPS(sites,L + (headoverride > 0 ? headoverride : QE))

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
  headoverride = para["headoverride"]

  if typeof(L) != Int
    L = L[1] * L[2]
  end 

  # set however many extra sites regarding the condition
  # experimental feature, QN AND QE
  # if we have both, then need AUX sites

  extras = headoverride > 0 ? 2 * headoverride :  QE * (QN + 1)

  #sites = Vector{Index}(undef, L + extras)
  sites = siteinds("Fermion", L + extras; conserve_qns =QN)
  # determines the number of total sites. IF QE and QN, then we have 2 sites for each QE, else 1
  # for s = 1: L + extras
  #   sites[s] =  siteind("Fermion"; addtags="n=$s", conserve_qns =QN)
  # end 

  #println(length(sites))
  return sites
end 


"""Preparing the state for the TE, if needed"""
function TE_stateprep(ψ, paras, sites)

  QE = paras["QE"]
  QN = paras["QN"]
  L = paras["L"]

  if QE == 0
    println("No QE, initial state as read")

  elseif QE == 1

    # if no QN, we add one 'ON' site as the left QE
    if !QN
      state = ["Occ"]

    # if QN, we attach two sites to the chain, with the QE site initiated as 'ON'
    else
      state = ["Emp", "Occ"]

    end 

    state = randomMPS( sites[1:QN + 1], state)
    ψ[1: QN + 1] = state

  elseif QE == 2
    
    if !QN
      state = ["Occ"]

    else
      state = ["Occ", "Emp"]
    end 

    state = randomMPS( sites[L + QN + 2: L + QN + 3 ], state)
    ψ[L + QN + 2: L+ QN + 3] = state

  else
    throw(ArgumentError("invalid QE number"))
  end 


  for i in eachindex(ψ)
    println(typeof(ψ[i]))
  end 

  return ψ

end 


function TE_stateprep(para)

  sites = init_site(para)
  L = para["L"]
  N = para["N"]
  QE = para["QE"]
  QN = para["QN"]

  if QN
    state = append!([ "Occ" for n=1:N] , ["Emp" for n=1:L -N]) 

    AUX1 = QE > 0 ? ["Emp"] : []
    QE1 = QE > 0 ? [ "Occ"] : []
    QE2 = QE > 1 ? [ "Occ"] : []
    AUX2 = QE > 1 ? ["Emp"] : []

    state = vcat(QE1, state)
    state = vcat(AUX1, state)
    state = vcat(state, QE2)
    state = vcat(state, AUX2)

    ψ = randomMPS(sites,state)

  else 
    ψ = randomMPS(sites,L + QE)
  end 

  
  return ψ, sites
end 