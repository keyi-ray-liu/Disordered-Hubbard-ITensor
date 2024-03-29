
"""add QE terms"""
function add_qe!(res, para::Dict,  L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1, which=1)

  println("Adding QE", which)
  QEen = para["QEen"][which]
  dp = para["dp"]
  ζ_dp = para["ζ_dp"]
  CN = para["CN"]
  QN = para["QN"]
  QEloc = para["QEloc"]
  scales = para["scales"]
  type = para["type"]
  scr = para["screening_qe"]
  range = para["range_qe"]

  Ltotal = get_systotal(para)

  # we had the check of dp length with QE
  dp = dp[which]
  ζ_dp = ζ_dp[which]


  # paux comes into play when QN is enabled, which is also encoded by head
  if which == 1
    paux = 1
    pqe = head 

  # for the rest of the QE, we currently assume they are all shuffled to the right, 
  else

    # aux sites are always to the right of QE for the RHS, i.e. QE2 at h + L + h, QE3 at h + L + h + h, etc.
    paux = Ltotal + head * which 
    pqe = Ltotal + head * which - 1
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

  else

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

function add_damping!(res, para, L, sites; head=0, which=2, QEmode="Damping-hop")

  println("Adding QE damping", which)
  QEen = para["QEen"][which]
  #QN = para["QN"]
  #QEloc = para["QEloc"]
  #scales = para["scales"]
  type = para["type"]
  #scr = para["screening_qe"]
  #range = para["range_qe"]

  Ltotal = get_systotal(para)


  if type == "Electron"
    opn = "Ntot"
    opc = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]

  elseif type == "Fermion"
    opn = "N"
    opc = [ ["C", "Cdag"]]
    
  end

  pexcur = head + Ltotal + ( which - 2) * 2 + 1
  pgscur = pexcur + 1
  pexnxt = pexcur + 2
  pgsnxt = pexcur + 3

  # we first bias the excited state. The add_qe function already biases QE2, so we only bias the nxt sites

  res += QEen, opn, pexnxt

  # we add specific interactions w.r.t the additional QEs

  if QEmode == "Damping-hop"
    for operator in opc

      opc1, opc2 = operator

      res += opc1, pexcur, opc2, pgscur, opc2, pexnxt, opc1, pgsnxt

      res += opc2, pgsnxt, opc1, pexnxt, opc1, pgscur, opc2, pexcur

    end

  elseif QEmode == "Damping-den"

    res += opn, pexcur, opn, pexnxt


  elseif QEmode == "Damping-reservoir"

    #overriding pcur because everything is coupled to QE2
    pexcur = head + Ltotal + 1
    pgscur = head + Ltotal + 2

    for operator in opc

      opc1, opc2 = operator


      res += which, opc1, pexcur, opc2, pgscur, opc2, pexnxt, opc1, pgsnxt

      res += which, opc2, pgsnxt, opc1, pexnxt, opc1, pgscur, opc2, pexcur

    end


  else

    error("Unrecognized damping mode!")
  end 

  return res
end 