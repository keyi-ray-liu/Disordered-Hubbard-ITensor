"""bulk hopping term"""
function add_hopping_bulk!(res, para::Dict, L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  println("Adding bulk hopping")

  decay = para["decay"]
  Ltotal = get_systotal(para)
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

  println("Adding SD hopping")
  t = para["t"]
  decay = para["decay"]
  Ltotal = get_systotal(para)
  type = para["type"]
  scales = para["scales"]
  
  sd_hop = para["sd_hop"]

  println(sd_hop)
  s_len = para["s_len"]
  d_len = para["d_len"]
  source_site = sd_hop["source_site"]
  drain_site = sd_hop["drain_site"]
  internal_hop = sd_hop["internal_hop"]
  to_chain_hop = sd_hop["to_chain_hop"]
  sd_loc = sd_hop["sd_loc"]
  sdcontact = sd_hop["sdcontact"]

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

  #determines if sd has contact
  if sdcontact
    p1 = head 
    p2 = Ltotal + head + 1
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

  if to_chain_hop != 0
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


  end 

  return res
end 