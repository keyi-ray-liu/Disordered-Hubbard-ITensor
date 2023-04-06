"""Hopping term"""
function add_hopping!(res, para::Dict, L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  t = para["t"]
  decay = para["decay"]
  Ltotal = prod(L)
  scales = para["scales"]

  # iterate through all sites to get NN
  for j=1  : Ltotal

    nns = get_nn(j, L)

    
    for nn in nns
      
      r = dis(j, nn, L, scales, disx, disy)
      hop = hopping(decay, r)

      #println( "$j, $nn, $hop")
      
      # Hopping
      p1 = j + head
      p2 = nn + head 
      s1 = sites[p1]
      s2 = sites[p2]

      if !if_gate
        res += -t * hop, "C", p1 ,"Cdag",p2
        res += -t * hop, "C", p2 ,"Cdag",p1

      else 
        hj =
        - t * hop * op("C", s1) * op("Cdag", s2) +
        - t * hop * op("C", s2) * op("Cdag", s1)

        gatefy!(res, factor, hj, τ)

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

  Ltotal = prod(L)

  for j=1 :Ltotal
    # E-E
    p1 = j + head
    s1 = sites[p1]

    for k= max(1, j - range):j-1

      # delta function setting up the exchange between nearest neighbor
      p2 = k + head
      s2 = sites[p2]

      # because k < j, we check if j is in the nn of k

      ifexch = (1 -  (j in get_nn(k, L)) * exch )
      r = dis(j, k, L, scales, disx, disy)
      
      #println("EE interaction: $j, $k, $ifexch, $r")
      if !if_gate
      # add the ee interaction term one by one
        res += λ_ee * ifexch / (  r + ζ_ee),"N",p1 ,"N", p2

      else
        eejk  = λ_ee * ifexch / ( r + ζ_ee) * op("N", s1) * op("N", s2)
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

  λ_ne = λ_ne * CN / Ltotal

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
      res += -cursum, "N", p1

    else 
      nej = -cursum * op("N", s1)
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

  if !if_gate
    res += QEen, "N", pqe

  else 
    g = QEen * op("N", sqe) 
    gatefy!(res, factor, g, τ)

  end 

  cavg = CN / Ltotal
  # dipole
  for i = 1 : Ltotal
    
    r = dis(i, QEloc[which], L, scales, disx, disy)

    r_dp = r ^ 3 + ζ_dp

    p1 = i + head
    s1 = sites[p1]


    if !if_gate
      if !QN
        res  +=  dp * r  / r_dp , "x", pqe, "N", p1
        # centering
        res  += -dp * r *  cavg / r_dp,  "x", pqe

      else
        res  +=  dp * r / r_dp, "C", paux, "Cdag", pqe, "N", p1
        res  +=  dp * r / r_dp, "C", pqe, "Cdag", paux, "N", p1
        
        # centering
        res  +=   - dp * r * cavg / r_dp, "C", paux, "Cdag", pqe
        res  +=   - dp * r * cavg / r_dp, "C", pqe, "Cdag", paux

      end 

    else

      if !QN
        dpla = dp * r  / r_dp * op( "x", sqe) * op("N", s1)
        dpl = -dp * r *  cavg / r_dp * op( "x", sqe)

      else

        #two level with AUX site
        dpla = dp * r / r_dp  * op( "C", saux) * op("Cdag", sqe) * op("N", s1) + 
                dp * r / r_dp  * op( "C", sqe) * op("Cdag", saux) * op("N", s1)

        # the offset term, to set the 'center of mass'
        # calculated as a uniform distribution:  L * N / 2
        dpl =  - dp * r * cavg / r_dp *  op("C", saux) * op( "Cdag", sqe)  
                - dp * r * cavg / r_dp *  op("C", sqe) * op( "Cdag", saux) 

      end 

      gatefy!(res, factor, dpla, τ)
      gatefy!(res, factor, dpl, τ)
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


