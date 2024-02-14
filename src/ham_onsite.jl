
"""add on-site interaction in the case of electrons"""
function add_onsite_hubbard!(res, para::Dict,  L::Vector{Int}, sites; if_gate=false, head=0, factor=2, τ=0.1)
  println("Adding on-site Hubbard")

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

  println("Adding on-site bias")
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
      res += bulk_bias[i], ops, p

    else 
      ons = bulk_bias[i] * op(ops, si)
      gatefy!(res, factor, ons, τ)
    end 
  end 

  return res
end 

"""add source-drain potential"""
function add_sd_potential(res, para, sites; if_gate = false, factor=2, τ=0.1)

  println("Adding SD bias")
  s_len = para["s_len"]
  d_len = para["d_len"]
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


