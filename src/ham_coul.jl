
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
    
    
    for k= max(1, j - range):j-1
      
      nns = allnn[k]
      # delta function setting up the exchange between nearest neighbor
      p2 = k + head
      s2 = sites[p2]

      # because k < j, we check if j is in the nn of k

      ifexch = (1 -  (j in keys(nns)) * exch )

      r = dis(j, k, L, scales, disx, disy)
      
      screen_factor = exp( - scr * r )
      
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


