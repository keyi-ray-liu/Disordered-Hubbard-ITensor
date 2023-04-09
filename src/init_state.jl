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
  mode = para["dynamode"]
  type = para["type"]

  
  # if not guess, using random MPS as init state
  if !guess
    Ltotal = prod(L)

    # Create an initial random matrix product state

    # QN and QE conditions if QN, explicit initiating conditions
    if QN

      # we already have check for the type of N
      if type == "Fermion"
        state = append!([ "Occ" for n=1:N] , ["Emp" for n=1:Ltotal -N]) 

      elseif type == "Electron"
        state = append!( ["Up" for n=1:N[1]], ["Dn" for n=1:N[2]], ["UpDn" for n=1:N[3]], ["Emp" for n=1: Ltotal - sum(N)])
      end 

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

        if mode == "both"
          front =  ["Emp", "Occ"] 
          back =  QE > 1 ? ["Occ", "Emp"] : []

        elseif mode == "left"
          front = ["Emp", "Occ"]
          back = QE > 1 ? front : []

        elseif mode == "right"
          front = ["Occ", "Emp"]
          back = QE > 1 ? front : []

        elseif mode == "empty"
          front  = ["Occ", "Emp"]
          back = QE > 1 ? ["Emp", "Occ"] : []

        elseif mode == "none"
          front =  ["Emp" for n=1:headoverride] 
          back = QE > 1 ? front : []

        else
          error("Unrecognized dyna mode")
        end 

        state = vcat(front, state)
        state = vcat(state, back)

      end 

      #@show state
      ψ0 = randomMPS(sites,state)

    # random MPS if QN is false, making sure no 'stuck' situation
    else

      ψ0 = randomMPS(sites, Ltotal + (headoverride > 0 ? headoverride : QE))

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
  type = para["type"]

  Ltotal = prod(L)

  # set however many extra sites regarding the condition
  # experimental feature, QN AND QE
  # if we have both, then need AUX sites

  extras = headoverride > 0 ? 2 * headoverride :  QE * (QN + 1)

  #sites = Vector{Index}(undef, L + extras)
  sites = siteinds(type, Ltotal + extras; conserve_qns =QN)
  # determines the number of total sites. IF QE and QN, then we have 2 sites for each QE, else 1
  # for s = 1: L + extras
  #   sites[s] =  siteind("Fermion"; addtags="n=$s", conserve_qns =QN)
  # end 

  #println(length(sites))
  return sites
end 


"""preparing a product state for TE"""
function TE_stateprep(para)

  sites = init_site(para)
  L = para["L"]
  N = para["N"]
  QE = para["QE"]
  QN = para["QN"]

  Ltotal = prod(L)

  if QN
    state = append!([ "Occ" for n=1:N] , ["Emp" for n=1: Ltotal -N]) 

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
    ψ = randomMPS(sites,Ltotal + QE)
  end 

  return ψ, sites
end 