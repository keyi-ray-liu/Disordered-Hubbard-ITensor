"""Generates the hamiltonian MPO for the given system configuration, both 1D and higher geometry"""
function init_ham(para::Dict,  L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; kwargs...)
  # parameters 

  println("H init start")

  factor= para["TEBDfactor"]
  τ = para["τ"]
  QE = para["QE"]
  QN = para["QN"]
  headoverride = para["headoverride"]
  type = para["type"]
  int_ee = para["int_ee"]
  int_ne = para["int_ne"]
  sd_override = para["sd_override"]
  U = para["U"]

  s_len = para["s_len"]
  d_len = para["d_len"]
  sd_hop = para["sd_hop"]


  if_gate = get(kwargs, :if_gate, false)
  
  mix_basis = get(sd_hop, "mix_basis", false)
  bulk_bias = get(sd_hop, "bulk_bias", 0)
  # if QE > 0, then at least left emitter, and if QN, we account for the AUX site

  # head denotes the position of QE1: 0 if QE == 0, 1 if QE but not QN, 2 if QE and QN
  # we add a head override function for the situation where we need 'empty' site indices
  head = headoverride > 0 ? headoverride : (QE > 0) * (QN + 1) + s_len
  println("head position is now at", head)

  if !if_gate
    res = OpSum()
  else
    res = ITensor[]
  end 

  
  ############################ begin chain hamiltonian ###################################

  
  res = add_hopping_bulk!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)

  # in case the 0 cases still adds overhead to hamiltonian
  if int_ee != 0

    res = add_ee!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)
  end 

  if int_ne != 0

    res = add_ne!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)
  end 

  # if we have spin in the system, add onsite interactions 
  if type == "Electron" && U != 0

    
    res = add_onsite_hubbard!(res, para, L, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)
  end 
  # ##################### end chain hamiltonian ###################################
  ##################### begin QE hamiltonian ###################################
  # QE part

  if headoverride == 0

    for which in 1:QE
      res = add_qe!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ, which=which)
    end
    
  end 


  ############################################## end QE hamiltonian ##########################################
  ######################### begin penalty ####################################

  # penalty terms for non-QN conserving hamiltonian
  # in the form of (sum_i n_i - N) ^2
  # if not QN, enforce
  if !QN 
    res = add_qn!(res, para, L, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)
  end 

  ########################### end penalty ###################################
  ##################### begin SD hamiltonian ###################################

  if s_len + d_len > 0 && sd_override == false

    if mix_basis
      
      energies = get(kwargs, :mix_basis_energies, [])
      ks = get(kwargs, :ks, [])
      LR = get(kwargs, :LR, [])

      if s_len + d_len != length(energies)
        error("Mix basis and SD(LR) does not match")
      end 

      if if_gate
        error("Does not support TEBD")
      end 
      # no TEBD support
      res = add_mix_sd(res, para, energies, ks, LR; head=head)


    else
      res = add_hopping_sd!(res, para, L, disx, disy, sites; if_gate = if_gate, head=head, factor =factor, τ=τ)
      res = add_sd_potential(res, para, sites; if_gate = if_gate,  factor=factor, τ=τ)
    end 
  end 

  if bulk_bias != 0
    res = add_onsite_bias!(res, para, sites, bulk_bias, if_gate=if_gate, head=head, factor= factor, τ=τ)
  end

  ########################### end SD hamiltonian ###################################

  #@show res 

  if !if_gate
    H = MPO(res, sites)

    @show maxlinkdim(H)
    return H

  else
    # reverse gates
    if factor == 2
      append!(res, reverse(res))
    end 
    
    return res

  end 
end


