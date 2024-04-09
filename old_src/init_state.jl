"""Returns a initial state.
Can either generate a random MPS or a guess free fermion initial state"""
function init_state(para, sites, disx, disy; kwargs...)

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
  s_len = para["s_len"]
  d_len = para["d_len"]
  
  source_config = get(kwargs, :source_config, [])
  drain_config = get(kwargs, :drain_config, [])
  state_override = get(kwargs, :state_override, [])

  if s_len != length(source_config) || d_len != length(drain_config)
    error("sd configs does not match inputs")
  end 

  Ltotal = get_systotal(para)
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

      if !isempty(state_override)

        
        state = [ op_str[i] for i in state_override]
        println("overriding init state: ", state)

      else
        state = [ op_str[i] for i= 1:4 for _ in 1:N[i] ]

      end

      if s_len + d_len > 0


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

        for _ in 2:QE
          state = vcat(state, QE2)
          state = vcat(state, AUX2)
        end 

      else

        if mode == "both"
          QEL =  [emp, occ ] 
          QER =  QE > 1 ? [occ, emp] : []

        elseif mode == "left"
          QEL = [emp, occ]
          QER = QE > 1 ? QEL : []

        elseif mode == "right"
          QEL = [occ, emp]
          QER = QE > 1 ? QEL : []

        elseif mode == "empty"
          QEL  = [occ, emp]
          QER = QE > 1 ? [emp, occ] : []

        elseif mode == "none"
          QEL =  [emp for n=1:headoverride] 
          QER = QE > 1 ? QEL : []

        else
          error("Unrecognized dyna mode")
        end 

        state = vcat(QEL, state)

        for _ in 2:QE
          state = vcat(state, QER)
        end

      end 


      @show length(state), length(sites)
      ψ0 = randomMPS(sites,state 
      #;linkdims=10
      )

    # random MPS if QN is false, making sure no 'stuck' situation
    else

      ψ0 = randomMPS(sites, Ltotal + (headoverride > 0 ? headoverride : QE)
      #; linkdims=10
      )

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


"""with QE X lattice"""
function init_state(para, sites, disx, disy, geometry::String; kwargs...)

  if geometry != "x-siam"

    return init_state(para, sites, disx, disy; kwargs...)

  end 

  L = para["L"]
  N = para["N"]
  QE = para["QE"]
  headoverride = para["headoverride"]
  mode = para["dynamode"]
  type = para["type"]
  state_override = get(kwargs, :state_override, [])

  half, leg = L

  #Ltotal = get_systotal(para)
  
  op_str = get_type_dict(type)
  emp = op_str[1]
  occ = op_str[2]
  
  # Create an initial random matrix product state

  # QN and QE conditions if QN, explicit initiating conditions


  # we already have check for the type of N

  if !isempty(state_override)

    
    state = [ op_str[i] for i in state_override]
    println("overriding init state: ", state)

  else
    state = [ op_str[i] for i= 1:4 for _ in 1:N[i] ]

  end

  Lstart = 1
  Lstate = [ state[Lstart + (h - 1) * leg: Lstart + h * leg - 1] for h in 1:half]

  midstart = 1 + half * leg
  midstate = [state[ midstart]]

  Rstart = midstart + 1
  Rstate = [ state[Rstart + (h -1) * leg: Rstart + h *leg - 1] for h in 1:half]


  if headoverride == 0

    QEL = QE > 0 ? [emp, occ] : []
    QER = QE >  half ? [ occ, emp ] : []

  else

    if mode == "both"
      QEL =  [emp, occ ] 
      QER =  QE > half ? [occ, emp] : []

    elseif mode == "left"
      QEL = [emp, occ]
      QER = QE > half ? QEL : []

    elseif mode == "right"
      QEL = [occ, emp]
      QER = QE > half ? QEL : []

    elseif mode == "empty"
      QEL  = [occ, emp]
      QER = QE > half ? [emp, occ] : []

    elseif mode == "none"
      QEL =  [emp for n=1:headoverride] 
      QER = QE > half ? QEL : []

    else
      error("Unrecognized dyna mode")
    end 

  end 

  # Lstates
  for h in 1:half
    state = vcat(state, QEL, Lstate[h])
  end 

  #mid
  state = vcat(state, midstate)

  #Rstate
  for h in 1:half
    state = vcat(state, Rstate[h], QER)
  end 


  @show length(state), length(sites)
  ψ0 = randomMPS(sites,state 
  #;linkdims=10
  )

# random MPS if QN is false, making sure no 'stuck' situation


  return ψ0
end

"""return the sites indices for further use, can add quantum emitters """
function init_site(para::Dict; kwargs...)
  L = para["L"]
  QE = para["QE"]
  QN = para["QN"]
  headoverride = para["headoverride"]
  type = para["type"]
  s_len = para["s_len"]
  d_len = para["d_len"]

  Ltotal = get_systotal(para)

  addtags = get(kwargs, :addtags, "")
  # set however many extra sites regarding the condition
  # experimental feature, QN AND QE
  # if we have both, then need AUX sites

  # we have previously check condition that SD and QE cannot both be greater than 0
  extras = headoverride > 0 ? QE * headoverride :  QE * (QN + 1) + s_len + d_len

  #sites = Vector{Index}(undef, L + extras)


  sites = siteinds(type, Ltotal + extras; conserve_qns =QN, addtags = addtags)

    #sites = siteinds( n -> n > s_len && n <= Ltotal + s_len ? type : "Boson", Ltotal + extras ; conserve_qns = QN )

  # determines the number of total sites. IF QE and QN, then we have 2 sites for each QE, else 1
  # for s = 1: L + extras
  #   sites[s] =  siteind("Fermion"; addtags="n=$s", conserve_qns =QN)
  # end 

  #println(length(sites))
  @show length(sites)
  return sites
end 


"""preparing a product state for TE"""
function TE_stateprep(para; kwargs...)

  sites = init_site(para; kwargs...)
  L = para["L"]
  N = para["N"]
  QE = para["QE"]
  QN = para["QN"]

  source_config = get(kwargs, :source_config, [])
  drain_config = get(kwargs, :drain_config, [])
  state_override = get(kwargs, :state_override, [])
  
  mode = para["dynamode"]
  Ltotal = get_systotal(para)
  type = para["type"]
  
  op_str = get_type_dict(type)
  emp = op_str[1]
  occ = op_str[2]

  s_len = para["s_len"]
  d_len = para["d_len"]

  if s_len != length(source_config) || d_len != length(drain_config)
    error("sd configs does not match inputs")
  end 

  
  if QN

    if !isempty(state_override)

      state = [ op_str[i] for i in state_override]

    else
      state = [ op_str[i] for i= 1:4 for _ in 1:N[i] ]

    end

    if s_len + d_len > 0

      source = [ op_str[n] for n in source_config]
      drain = [ op_str[n] for n in drain_config]
      state = vcat(source, state, drain)

    else

      if mode == "both"
        QEL =  [emp, occ ] 
        QER =  QE > 1 ? [occ, emp] : []

      elseif mode == "left"
        QEL = [emp, occ]
        QER = QE > 1 ? QEL : []

      elseif mode == "right"
        QEL = [occ, emp]
        QER = QE > 1 ? QEL : []

      elseif mode == "empty"
        QEL  = [occ, emp]
        QER = QE > 1 ? [emp, occ] : []

      elseif mode == "none"
        QEL =  [emp for n=1:headoverride] 
        QER = QE > 1 ? QEL : []

      else
        error("Unrecognized dyna mode")
      end 

      state = vcat(QEL, state)
      state = vcat(state, QER)
    end

    print(state)
    ψ = randomMPS(sites,state
    #; linkdims=10
    )

  else 
    ψ = randomMPS(sites,Ltotal + QE
    #; linkdims=10
    )
  end 

  return ψ, sites
end 