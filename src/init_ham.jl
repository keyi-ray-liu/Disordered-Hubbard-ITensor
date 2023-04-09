"""Generates the hamiltonian MPO for the given system configuration, both 1D and higher geometry"""
function init_ham(para::Dict,  L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false)
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
  
  # if QE > 0, then at least left emitter, and if QN, we account for the AUX site

  # head denotes the position of QE1: 0 if QE == 0, 1 if QE but not QN, 2 if QE and QN
  # we add a head override function for the situation where we need 'empty' site indices
  head = headoverride > 0 ? headoverride : (QE > 0) * (QN + 1) 
  println("head position is now at", head)

  if !if_gate
    res = OpSum()
  else
    res = ITensor[]
  end 

  
  ############################ begin chain hamiltonian ###################################
  res = add_hopping!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)

  # in case the 0 cases still adds overhead to hamiltonian
  if int_ee != 0
    res = add_ee!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)
  end 

  if int_ne != 0
    res = add_ne!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)
  end 

  # if we have spin in the system, add onsite interactions 
  if type == "Electron"
    res = add_onsite!(res, para, L, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)
  end 
  # ##################### end chain hamiltonian ###################################
  ##################### begin QE hamiltonian ###################################
  # QE part
  # left QE
  if QE > 0 && headoverride == 0
    res = add_qe!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ, which=1)
  end 

  # right QE
  if QE > 1 && headoverride == 0
    res = add_qe!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ, which=2)
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
  if !if_gate
    H = MPO(res, sites)
    return H

  else
    # reverse gates
    if factor == 2
      append!(res, reverse(res))
    end 
    
    return res

  end 
end


