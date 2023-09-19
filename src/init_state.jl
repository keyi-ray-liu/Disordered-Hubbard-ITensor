"""Returns a initial state.
Can either generate a random MPS or a guess free fermion initial state"""
function init_state(para, sites, disx, disy)

  L = para["L"]
  N = para["N"]
  guess = para["guess"]
  t = para["t"]
  λ_ne = para["int_ne"]
  ζ_ne, _ = para["ζ"]
  decay = para["decay"]
  self = para["self_nuc"]
  QE = para["QE"]
  QN = para["QN"]
  headoverride = para["headoverride"]
  mode = para["dynamode"]
  type = para["type"]
  
  source_config = para["source_config"]
  drain_config = para["drain_config"]

  Ltotal = prod(L)
  scales = para["scales"]
  snake = para["snake"]
  QEen = para["QEen"]
  
  op_str = get_type_dict(type)
  emp = op_str[1]
  occ = op_str[2]

  # if not guess, using random MPS as init state
  if !guess
    

    # Create an initial random matrix product state

    # QN and QE conditions if QN, explicit initiating conditions
    if QN

      # we already have check for the type of N

      state = [ op_str[i] for i= 1:4 for _ in 1:N[i] ]

      if length(source_config) + length(drain_config) > 0

        source = [ op_str[n] for n in source_config]
        drain = [ op_str[n] for n in drain_config]

        state = vcat(source, state, drain)

      elseif headoverride == 0

        AUX1 = QE > 0 ? [emp] : []
        QE1 = QE > 0 ? [ occ ] : []
        QE2 = QE > 1 ? [ occ ] : []
        AUX2 = QE > 1 ? [emp] : []

        state = vcat(QE1, state)
        state = vcat(AUX1, state)
        state = vcat(state, QE2)
        state = vcat(state, AUX2)

      else

        if mode == "both"
          front =  [emp, occ ] 
          back =  QE > 1 ? [occ, emp] : []

        elseif mode == "left"
          front = [emp, occ]
          back = QE > 1 ? front : []

        elseif mode == "right"
          front = [occ, emp]
          back = QE > 1 ? front : []

        elseif mode == "empty"
          front  = [occ, emp]
          back = QE > 1 ? [emp, occ] : []

        elseif mode == "none"
          front =  [emp for n=1:headoverride] 
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
    

    # single particle H
    H = zeros((Ltotal + QE * (QN + 1), Ltotal + QE * (QN + 1)))

    #hopping
    for i = 1:Ltotal

      nns = get_nn(i, L, snake =snake)

      for nn in nns
        r = dis(i, nn, L, scales, disx, disy)
        hop = disorder_hopping(decay, r)
        
        H[i + (QN + 1), nn + (QN + 1)] = -t * hop
        H[nn + (QN + 1), i + (QN + 1)] = -t * hop
      end 

    end 
    
    # nuclear

    if λ_ne != 0
      for j = 1:Ltotal
        if self
          cursum = 0
        else
          cursum = -λ_ne / ζ_ne
        end

        for k=1:Ltotal
          cursum += λ_ne / ( dis(j, k, L, scales, disx, disy) + ζ_ne)
        end

        H[j + (QN + 1), j + (QN + 1)] = -cursum

      end 
    end 

    if QE > 0

      H[ QN + 1, QN + 1] = QEen

    end 

    if QE > 1

      H[ Ltotal  + QN + 2, Ltotal + QN + 2] = QEen

    end 

    @show ishermitian(H)
    _, u = eigen(H)

  
    # Get the Slater determinant
    # Create an mps for the free fermion ground state

    if type == "Fermion"
      ϕ = u[:, 1:sum(N)]
      ψ0 = slater_determinant_to_mps(sites, ϕ)

    elseif type == "Electron"

      @show N

      ϕup = u[:, 1:N[1]]
      ϕdn = u[:, 1:N[2]]

      ψ0 = slater_determinant_to_mps(sites, ϕup, ϕdn, eigval_cutoff=1E-10, maxblocksize=4)

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
  source_config = para["source_config"]
  drain_config = para["drain_config"]
  LSR_bruteforce = para["LSR_bruteforce"]

  Ltotal = prod(L)

  # set however many extra sites regarding the condition
  # experimental feature, QN AND QE
  # if we have both, then need AUX sites

  # we have previously check condition that SD and QE cannot both be greater than 0
  extras = headoverride > 0 ? 2 * headoverride :  QE * (QN + 1) + length(source_config) + length(drain_config)

  #sites = Vector{Index}(undef, L + extras)

  if !LSR_bruteforce
    sites = siteinds(type, Ltotal + extras; conserve_qns =QN)

  else
    sites = siteinds( n -> n > length(source_config) && n <= Ltotal + length(source_config) ? type : "Boson", Ltotal + extras ; conserve_qns = QN )
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

  source_config = para["source_config"]
  drain_config = para["drain_config"]

  Ltotal = prod(L)
  type = para["type"]
  
  op_str = get_type_dict(type)
  emp = op_str[1]
  occ = op_str[2]
  
  if QN
    state = [ op_str[i] for i= 1:4 for _ in 1:N[i] ]

    if length(source_config) + length(drain_config) > 0

      source = [ op_str[n] for n in source_config]
      drain = [ op_str[n] for n in drain_config]
      state = vcat(source, state, drain)

    else
      AUX1 = QE > 0 ? [emp] : []
      QE1 = QE > 0 ? [ occ] : []
      QE2 = QE > 1 ? [ occ] : []
      AUX2 = QE > 1 ? [emp] : []

      state = vcat(QE1, state)
      state = vcat(AUX1, state)
      state = vcat(state, QE2)
      state = vcat(state, AUX2)
    end

    print(state)
    ψ = randomMPS(sites,state)

  else 
    ψ = randomMPS(sites,Ltotal + QE)
  end 

  return ψ, sites
end 