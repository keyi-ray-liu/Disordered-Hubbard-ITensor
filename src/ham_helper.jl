"""bulk hopping term"""
function add_hopping_bulk!(res, para::Dict, L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  decay = para["decay"]
  Ltotal = prod(L)
  scales = para["scales"]
  type = para["type"]
  snake = para["snake"]
  allnn = para["allnn"]
  # iterate through all sites to get NN

  if type == "Fermion"
    operators = [ ["C", "Cdag"]]

  elseif type == "Electron"
    operators = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]

  end 

  for j=1  : Ltotal

    nns = allnn[j]
    for nn in keys(nns)

      t = nns[nn]
      
      r = dis(j, nn, L, scales, disx, disy)
      hop = disorder_hopping(decay, r) * t

      # temporary fix! 
      if snake
        hop = 1.0
      end 

      println( "adding hopping $(j + head), $(nn + head), $( hop)")
      
      # Hopping
      p1 = j + head
      p2 = nn + head 
      s1 = sites[p1]
      s2 = sites[p2]

      for operator in operators

        op1, op2 = operator 

        if !if_gate
          res += hop, op1, p1 ,op2, p2
          res += hop, op1, p2 ,op2, p1

        else 
          hj =
          hop * op(op1, s1) * op(op2, s2) +
          hop * op(op1, s2) * op(op2, s1)

          gatefy!(res, factor, hj, τ)

        end 

      end 

    end 

  end

  return res
end 


"""sd hopping term"""
function add_hopping_sd!(res, para::Dict, L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  t = para["t"]
  decay = para["decay"]
  Ltotal = prod(L)
  type = para["type"]
  scales = para["scales"]
  
  sd_hop = para["sd_hop"]

  println(sd_hop)
  s_len = length(para["source_config"])
  d_len = length(para["drain_config"])
  source_site = sd_hop["source_site"]
  drain_site = sd_hop["drain_site"]
  internal_hop = sd_hop["internal_hop"]
  to_chain_hop = sd_hop["to_chain_hop"]
  sd_loc = sd_hop["sd_loc"]

  # iterate through all sites to get NN

  if type == "Fermion"
    bulk_operators = [ ["C", "Cdag"]]

  elseif type == "Electron"
    bulk_operators = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]

  end 


  sd_operators = bulk_operators

  #source, internal
  for i_source = 1  : s_len - 1

    #println("internal hop, source, $i_source")
    #nns = get_nn(j, L, snake=snake, geometry=geometry, belong="source")
    # for the moment we only consider 1D SD

    p1 = i_source
    p2 = i_source + 1
    s1 = sites[p1]
    s2 = sites[p2]

    for operator in sd_operators

      op1, op2 = operator 

      if !if_gate
        res += t * internal_hop, op1, p1 ,op2, p2
        res += t * internal_hop, op1, p2 ,op2, p1

      else 
        hj =
        t  * internal_hop *  op(op1, s1) * op(op2, s2) +
        t  * internal_hop * op(op1, s2) * op(op2, s1)

        gatefy!(res, factor, hj, τ)

      end 

    end 

  end

  #drain, internal
  for i_drain =1  : d_len - 1

    #nns = get_nn(j, L, snake=snake, geometry=geometry, belong="source")
    # for the moment we only consider 1D SD

    #println("internal hop, drain, $i_drain")
    p1 = i_drain + head + Ltotal 
    p2 = i_drain + head + Ltotal + 1
    s1 = sites[p1]
    s2 = sites[p2]

    for operator in sd_operators

      op1, op2 = operator 

      if !if_gate
        res += t * internal_hop , op1, p1 ,op2, p2
        res += t * internal_hop, op1, p2 ,op2, p1

      else 
        hj =
        t  * internal_hop * op(op1, s1) * op(op2, s2) +
        t  * internal_hop * op(op1, s2) * op(op2, s1)

        gatefy!(res, factor, hj, τ)

      end 

    end 

  end

  #source, onto chain
  # 
  r_source = dis(source_site, sd_loc[1], L, scales, disx, disy)
  hop = disorder_hopping(decay, r_source) * to_chain_hop

  ps = head
  pb = source_site + head 
  ss = sites[ps]
  sb = sites[pb]

  
  for (k, sd_operator) in enumerate(sd_operators)
    
    sop1, sop2 = sd_operator
    bop1, bop2 = bulk_operators[k]

    #println("source hop to bulk, from $ps, $sop1, to $pb, $bop2, $$hop")
    #println("source hop to bulk, from $pb, $bop2, to $pb, $bop2, $$hop")
    if !if_gate
      res += t * hop, sop1, ps ,bop2, pb
      res += t * hop, bop1, pb ,sop2, ps

    else 
      hj =
      t * hop * op(sop1, ss) * op(bop2, sb) +
      t * hop * op(bop1, sb) * op(sop2, ss)

      gatefy!(res, factor, hj, τ)
    end 



  end 

  # drain, onto chain]

  r_drain = dis(drain_site, sd_loc[2], L, scales, disx, disy)
  hop = disorder_hopping(decay, r_drain) * to_chain_hop

  pb = drain_site + head
  pd = Ltotal + head + 1

  #println("drain hop to bulk, $p1, $p2, $hop")
  sb = sites[pb]
  sd = sites[pd]

  for (k, sd_operator) in enumerate(sd_operators)
    
    dop1, dop2 = sd_operator
    bop1, bop2 = bulk_operators[k]
    if !if_gate
        res += t * hop, dop1, pd ,bop2, pb
        res += t * hop, bop1, pb ,dop2, pd

    else 
      hj =
      t * hop * op(dop1, sd) * op(bop2, sb) +
      t * hop * op(bop1, sb) * op(dop2, sd)

      gatefy!(res, factor, hj, τ)
    end 


  end 

  return res
end 
"""E-E interaction term"""
function add_ee!(res, para::Dict,  L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  λ_ee = para["int_ee"]
  ζ_ee = para["ζ"][2]
  exch = para["exch"]
  range = para["range"]
  scales = para["scales"]
  #snake = para["snake"]
  scr = para["screening_int"]
  type = para["type"]
  Ltotal = prod(L)
  allnn = para["allnn"]

  if type == "Fermion"
    ops = "N"

  elseif type == "Electron"
    ops = "Ntot"

  end 

  for j=1 :Ltotal
    # E-E
    p1 = j + head
    s1 = sites[p1]
    nns = allnn[j]
    
    for k= max(1, j - range):j-1
      
      # delta function setting up the exchange between nearest neighbor
      p2 = k + head
      s2 = sites[p2]

      # because k < j, we check if j is in the nn of k

      ifexch = (1 -  (j in keys(nns)) * exch )
      r = dis(j, k, L, scales, disx, disy)
      
      screen_factor = exp( - scr * r )
      #println("EE interaction: $j, $k, $ifexch, $r")
      if !if_gate
      # add the ee interaction term one by one
        res += λ_ee * ifexch * screen_factor / (  r + ζ_ee), ops,p1 ,ops, p2

      else
        eejk  = λ_ee * ifexch * screen_factor / ( r + ζ_ee) * op(ops, s1) * op(ops, s2)
        gatefy!(res, factor, eejk, τ)

      end 

    end

  end 

  return res

end 

"""N-E interaction term"""
function add_ne!(res, para::Dict,  L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  λ_ne = para["int_ne"]
  ζ_ne = para["ζ"][1]
  self = para["self_nuc"]
  CN = para["CN"]
  Ltotal = prod(L)
  scales = para["scales"]
  scr = para["screening_int"]
  type = para["type"]
  range = para["range"]

  λ_ne = λ_ne * CN / Ltotal

  if type == "Fermion"
    ops = "N"

  elseif type == "Electron"
    ops = "Ntot"

  end 

  for j=1 :Ltotal

    p1 = j + head
    s1 = sites[p1]
    
    # N-E

    # check if self self_nuc

    if self
      cursum = 0
    else
      cursum = -λ_ne / ζ_ne
    end

    for l= max(1, j - range) : min(Ltotal, j + range)

      r = dis(j, l, L, scales, disx, disy)

      screen_factor = exp( - scr * r )
      #println("NE interaction: $j, $l, $r")
      cursum += λ_ne * screen_factor / ( r + ζ_ne)
    end

    if !if_gate
      res += -cursum, ops, p1

    else 
      nej = -cursum * op(ops, s1)
      gatefy!(res, factor, nej, τ)
    end 

  end

  return res

end 

"""add QE terms"""
function add_qe!(res, para::Dict,  L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1, which=1)

  QEen = para["QEen"]
  dp = para["dp"]
  ζ_dp = para["ζ_dp"]
  CN = para["CN"]
  QN = para["QN"]
  QEloc = para["QEloc"]
  scales = para["scales"]
  type = para["type"]
  scr = para["screening_qe"]
  range = para["range_qe"]

  Ltotal = prod(L)

  # we had the check of dp length with QE
  dp = dp[which]
  ζ_dp = ζ_dp[which]

  if which == 1
    paux = 1
    pqe = head 

  elseif which == 2

    paux = Ltotal + head * 2
    pqe = Ltotal + head + 1
  end 

  saux = sites[paux]
  sqe = sites[pqe]

  # diagonal energy term
  # for the spinless case, the occ state of QE pos takes the 'ex' state 
  # for the spinful case, the QE are initialized in 'up', and cdagup cup takes the charge between the ex and gs of the QE, the nhat is replaced by Ntotal

  # set operators based on the system
  if type == "Electron"
    opn = "Ntot"
    opc = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]

  elseif type == "Fermion"
    opn = "N"
    opc = [ ["C", "Cdag"]]
    
  end

  if !if_gate
    res += QEen, opn, pqe

  else 
    g = QEen * op(opn, sqe) 
    gatefy!(res, factor, g, τ)
  end 

  cavg = CN / Ltotal
  # dipole
  if which == 1

    qe_begin = 1
    qe_end = min(range + 1, Ltotal)

  elseif which == 2

    qe_begin = max(Ltotal - range , 1)
    qe_end = Ltotal 
  end 

  for i = qe_begin : qe_end
    
    r = dis(i, QEloc[which], L, scales, disx, disy)

    screen_factor = exp( - scr * r)
    println("distance to QE: $r")
    r_dp = r ^ 3 + ζ_dp

    p1 = i + head
    s1 = sites[p1]

    for operator in opc

      opc1, opc2 = operator

      if !if_gate
        if !QN
          res  +=  dp * r * screen_factor / r_dp , "x", pqe, opn, p1
          # centering
          res  += -dp * r * screen_factor * cavg / r_dp,  "x", pqe

        else
          res  +=  dp * r * screen_factor / r_dp, opc1, paux, opc2, pqe, opn, p1
          res  +=  dp * r * screen_factor / r_dp, opc1, pqe, opc2, paux, opn, p1
          
          # centering
          res  +=   - dp * r * screen_factor * cavg / r_dp, opc1, paux, opc2, pqe
          res  +=   - dp * r * screen_factor * cavg / r_dp, opc1, pqe, opc2, paux

        end 

      else

        if !QN
          dpla = dp * r * screen_factor / r_dp * op( "x", sqe) * op(opn, s1)
          dpl = -dp * r * screen_factor * cavg / r_dp * op( "x", sqe)

        else

          #two level with AUX site
          dpla = dp * r * screen_factor / r_dp  * op( opc1, saux) * op(opc2, sqe) * op(opn, s1) + 
                  dp * r * screen_factor / r_dp  * op( opc1, sqe) * op(opc2, saux) * op(opn, s1)

          # the offset term, to set the 'center of mass'
          # calculated as a uniform distribution:  L * N / 2
          dpl =  - dp * r * screen_factor * cavg / r_dp *  op(opc1, saux) * op( opc2, sqe)  
                  - dp * r * screen_factor * cavg / r_dp *  op(opc1, sqe) * op( opc2, saux) 

        end 

        gatefy!(res, factor, dpla, τ)
        gatefy!(res, factor, dpl, τ)
      end 

    end 


    # the offset term, to set the 'center of mass'
    # calculated as a uniform distribution:  L * N / 2
    #res +=  dp * L * N / ( 2 * r^3), "x", head, "N", left + head
  end 

  return res

end 

"""add QN terms if necessary"""
function add_qn!(res, para::Dict,  L::Vector{Int}, sites; if_gate=false, head=0, factor=2, τ=0.1, Λ=30)
  N = para["N"]
  Ltotal = prod(L)

  for i= 1:Ltotal

    p1 = i + head
    s1 = sites[p1]
    # linear terms

    if !if_gate
      res += - 2 * Λ * N, "N", s1 + head
    else
      li = - 2 * Λ * N * op( "N", s1)
      gatefy!(res, factor, li, τ)
    end 
      

    # quadratic terms
    for j =1:Ltotal

      p2 = j + head
      s2 = sites[p2]

      if !if_gate
        res += Λ, "N", p1, "N", p2
      else

        if s1 != s2
          lij = Λ * op("N", s1) * op("N", s2)
        end 
        gatefy!(res, factor, lij, τ)

      end 
    end 

  end 

  return res

end 

"""add on-site interaction in the case of electrons"""
function add_onsite_hubbard!(res, para::Dict,  L::Vector{Int}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  U = para["U"]
  Ltotal = prod(L)

  for i = 1:Ltotal

    p = i + head
    si = sites[p]

    if !if_gate
      res += U, "Nupdn", p

    else 
      ons = U * op("Nupdn", si)
      gatefy!(res, factor, ons, τ)
    end 
  end 

  return res
end 


""" adding bulk bias"""
function add_onsite_bias!(res, para::Dict, sites, bulk_bias; if_gate=false, head=0, factor=2, τ=0.1)

  L = para["L"]
  Ltotal = prod(L)
  type = para["type"]

  if type == "Fermion"
    ops = "N"

  elseif type == "Electron"
    ops = "Ntot"
  end 

  for i = 1:Ltotal

    #println("adding bulk bias, $(i + head), $bulk_bias")
    p = i + head
    si = sites[p]

    if !if_gate
      res += bulk_bias, ops, p

    else 
      ons = bulk_bias * op(ops, si)
      gatefy!(res, factor, ons, τ)
    end 
  end 

  return res
end 

"""add source-drain potential"""
function add_sd_potential(res, para, sites; if_gate = false, factor=2, τ=0.1)

  s_len = length(para["source_config"])
  d_len = length(para["drain_config"])
  sd_hop = para["sd_hop"]
  L = para["L"]
  type = para["type"]
  Ltotal = prod(L)

  if type == "Fermion" 
    ops = "N"

  elseif type == "Electron"
    ops = "Ntot"
  end 

  source_offset = sd_hop["source_offset"]
  drain_offset = sd_hop["drain_offset"]

  for i_source in 1:s_len

    ps = i_source
    ss = sites[ps]

    #println("adding source offset, $ps, $source_offset, $ops")
    if !if_gate
      res += source_offset, ops, ps

    else 
      ons = source_offset * op(ops, ss)
      gatefy!(res, factor, ons, τ)
    end 

  end 

  for i_drain in 1:d_len

    pd = i_drain + s_len + Ltotal
    sd = sites[pd]

    #println("adding drain offset, $pd, $drain_offset. $ops")

    if !if_gate
      res += drain_offset, ops, pd

    else 
      ons = drain_offset * op(ops, sd)
      gatefy!(res, factor, ons, τ)
    end 

  end 

  return res

end 


"""add mixed basis hopping to chain"""
function add_mix_sd(res, para, energies, ks; head=0)


  sd_hop = para["sd_hop"]
  t = para["t"]
  s_len = length(para["source_config"])
  d_len = length(para["drain_config"])
  L = para["L"]
  type = para["type"]

  source_site = sd_hop["source_site"]
  drain_site = sd_hop["drain_site"]

  to_chain_hop = sd_hop["to_chain_hop"]

  Ltotal = prod(L)

  if type == "Fermion"
    ops = [ ["Cdag", "C"]]
    opn = "N"

  elseif type == "Electron"
    ops = [ ["Cdagup", "Cup"], ["Cdagdn", "Cdn"]]
    opn = "Ntot"

  end 
  # source, i.e. L
  for k in 1:s_len

    mixed_k = ks[k]
    hopping = t * to_chain_hop * Ukj(mixed_k, 1, s_len)

    for (op1, op2) in ops
      # hop to system
      res +=  hopping, op1, head + source_site, op2, k
      res +=  hopping, op1, k, op2, head + source_site
    end

      # diagonal
    res += energies[k], opn, k
  end 

  # drain, i.e. R
  for k in 1:d_len

    mixed_k = ks[k + s_len]
    hopping = t * to_chain_hop * Ukj(mixed_k, 1, d_len)

    for (op1, op2) in ops

      # hop to system
      res +=  hopping, op1, head + drain_site, op2, k + head + Ltotal
      res +=  hopping, op1, k + head + Ltotal, op2, head + drain_site

    end 

    # diagonal
    res += energies[k + s_len], opn, k + head + Ltotal
    
  end 

  return res

  
end

