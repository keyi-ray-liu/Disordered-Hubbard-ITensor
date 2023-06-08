"""Hopping term"""
function add_hopping!(res, para::Dict, L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  t = para["t"]
  decay = para["decay"]
  Ltotal = prod(L)
  scales = para["scales"]
  type = para["type"]
  snake = para["snake"]
  geometry = para["geometry"]
  
  # iterate through all sites to get NN

  if type == "Fermion"
    operators = [ ["C", "Cdag"]]

  elseif type == "Electron"
    operators = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]

  end 

  for j=1  : Ltotal

    nns = get_nn(j, L, snake=snake, geometry=geometry)

    println("site $j, NN $nns")

    for nn in nns
      
      r = dis(j, nn, L, scales, disx, disy)
      hop = hopping(decay, r)

      # temporary fix! 
      if snake
        hop = 1.0
      end 

      #println( "$j, $nn, $hop")
      
      # Hopping
      p1 = j + head
      p2 = nn + head 
      s1 = sites[p1]
      s2 = sites[p2]

      for operator in operators

        op1, op2 = operator 

        if !if_gate
          res += -t * hop, op1, p1 ,op2, p2
          res += -t * hop, op1, p2 ,op2, p1

        else 
          hj =
          - t * hop * op(op1, s1) * op(op2, s2) +
          - t * hop * op(op1, s2) * op(op2, s1)

          gatefy!(res, factor, hj, τ)

        end 

      end 

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
  snake = para["snake"]

  type = para["type"]
  Ltotal = prod(L)

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

      # delta function setting up the exchange between nearest neighbor
      p2 = k + head
      s2 = sites[p2]

      # because k < j, we check if j is in the nn of k

      ifexch = (1 -  (j in get_nn(k, L, snake=snake)) * exch )
      r = dis(j, k, L, scales, disx, disy)
      
      #println("EE interaction: $j, $k, $ifexch, $r")
      if !if_gate
      # add the ee interaction term one by one
        res += λ_ee * ifexch / (  r + ζ_ee), ops,p1 ,ops, p2

      else
        eejk  = λ_ee * ifexch / ( r + ζ_ee) * op(ops, s1) * op(ops, s2)
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

  type = para["type"]

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

    for l=1 :Ltotal

      r = dis(j, l, L, scales, disx, disy)
      #println("NE interaction: $j, $l, $r")
      cursum += λ_ne / ( r + ζ_ne)
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
  for i = 1 : Ltotal
    
    r = dis(i, QEloc[which], L, scales, disx, disy)

    println("distance to QE: $r")
    r_dp = r ^ 3 + ζ_dp

    p1 = i + head
    s1 = sites[p1]

    for operator in opc

      opc1, opc2 = operator

      if !if_gate
        if !QN
          res  +=  dp * r  / r_dp , "x", pqe, opn, p1
          # centering
          res  += -dp * r *  cavg / r_dp,  "x", pqe

        else
          res  +=  dp * r / r_dp, opc1, paux, opc2, pqe, opn, p1
          res  +=  dp * r / r_dp, opc1, pqe, opc2, paux, opn, p1
          
          # centering
          res  +=   - dp * r * cavg / r_dp, opc1, paux, opc2, pqe
          res  +=   - dp * r * cavg / r_dp, opc1, pqe, opc2, paux

        end 

      else

        if !QN
          dpla = dp * r  / r_dp * op( "x", sqe) * op(opn, s1)
          dpl = -dp * r *  cavg / r_dp * op( "x", sqe)

        else

          #two level with AUX site
          dpla = dp * r / r_dp  * op( opc1, saux) * op(opc2, sqe) * op(opn, s1) + 
                  dp * r / r_dp  * op( opc1, sqe) * op(opc2, saux) * op(opn, s1)

          # the offset term, to set the 'center of mass'
          # calculated as a uniform distribution:  L * N / 2
          dpl =  - dp * r * cavg / r_dp *  op(opc1, saux) * op( opc2, sqe)  
                  - dp * r * cavg / r_dp *  op(opc1, sqe) * op( opc2, saux) 

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
function add_onsite!(res, para::Dict,  L::Vector{Int}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  U = para["U"]
  Ltotal = prod(L)

  for i = 1:Ltotal

    pi = i + head
    si = sites[pi]

    if !if_gate
      res += U, "Nupdn", pi

    else 
      ons = U * op("Nupdn", si)
      gatefy!(res, factor, ons, τ)
    end 
  end 

  return res
end 

