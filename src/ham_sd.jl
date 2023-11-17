
"""add mixed basis hopping to chain"""
function add_mix_sd(res, para, energies, ks, LR; head=0)

  offsetLR(v) = v > s_len ? Ltotal : 0
  sd_hop = para["sd_hop"]
  t = para["t"]
  s_len = para["s_len"]
  d_len = para["d_len"]
  L = para["L"]
  type = para["type"]

  source_site = sd_hop["source_site"]
  drain_site = sd_hop["drain_site"]

  to_chain_hop = sd_hop["to_chain_hop"]
  interacting = sd_hop["interacting"]
  inelastic = sd_hop["inelastic"]
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

  # testing inelastic, we add a hopping between the nearest neighbor mode between L and r

  if inelastic
    for k1 in 1: (s_len + d_len)

      key = LR[k1] == 1 ? 2 : 1
      diff = findfirst(==(key), LR[ k1 + 1:end])

      
      # we add a pair of interaction term
      if !isnothing(diff)
        
        k2 = k1 + diff

        for (op1, op2) in ops
          hopping = t  * exp( - 2.0 * abs( energies[k1] - energies[k2]))
          
          res += hopping, op1, k1 + offsetLR(k1), op2, k2 + offsetLR(k2)
          res += hopping, op1, k2 + offsetLR(k2) , op2, k1 + offsetLR(k1)
        end
      end 

    end 
  end

  return res

  
end

