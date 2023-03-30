"""Hopping term"""
function add_hopping!(res, para::Dict, L::Int, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  t = para["t"]
  decay = para["decay"]

  # adjust the site based on if there are left emitter
  for j=1  :L-1 

    r = dis(j, j + 1, disx, disy)
    hop = hopping(decay, r)

    #println(hop)
    # Hopping
    p1 = j + head
    p2 = j + head + 1
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

  return res
end 

"""E-E interaction term"""
function add_ee!(res, para::Dict, L::Int, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  λ_ee = para["int_ee"]
  ζ_ee = para["ζ"][2]
  exch = para["exch"]
  range = para["range"]

  for j=1 :L 
    # E-E
    p1 = j + head
    s1 = sites[p1]

    for k= max(1, j - range):j-1

      # delta function setting up the exchange between nearest neighbor
      p2 = k + head
      s2 = sites[p2]

      ifexch = (1 -  ==(1, abs(j - k)) * exch )
      
      if !if_gate
      # add the ee interaction term one by one
        res += λ_ee * ifexch / ( dis(j, k, disx, disy) + ζ_ee),"N",p1 ,"N", p2

      else
        eejk  = λ_ee * ifexch / ( dis(j, k, disx, disy) + ζ_ee) * op("N", s1) * op("N", s2)
        gatefy!(res, factor, eejk, τ)

      end 

    end

  end 

  return res

end 


"""N-E interaction term"""
function add_ne!(res, para::Dict, L::Int, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1)

  λ_ne = para["int_ne"]
  ζ_ne = para["ζ"][1]
  self = para["self_nuc"]
  CN = para["CN"]

  λ_ne = λ_ne * CN / L

  for j=1 :L 

    p1 = j + head
    s1 = sites[p1]
    
    # N-E

    # check if self self_nuc

    if self
      cursum = 0
    else
      cursum = -λ_ne / ζ_ne
    end

    for l=1 :L 
      cursum += λ_ne / ( dis(j, l, disx, disy) + ζ_ne)
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
function add_qe!(res, para::Dict, L::Int, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false, head=0, factor=2, τ=0.1, which=1)

  QEen = para["QEen"]
  dp = para["dp"]
  ζ_dp = para["ζ_dp"]
  QEoffset = para["QEoffset"]
  CN = para["CN"]
  QN = para["QN"]
  QE = para["QE"]


  # we had the check of dp length with QE
  dp = dp[which]
  ζ_dp = ζ_dp[which]

  if which == 1
    paux = 1
    pqe = head 

  elseif which == 2

    paux = L + head * 2
    pqe = L + head + 1
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

  cavg = CN / L 
  # dipole
  for i = 1 : L 
    
    if which == 1
      r = dis(i, QEoffset, disx, disy)

    elseif which == 2
      r = L - 1 + QE * (QEoffset + 1) - dis(i, QEoffset, disx, disy)
    end 

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
function add_qn!(res, para::Dict, L::Int, sites; if_gate=false, head=0, factor=2, τ=0.1, Λ=30)
  N = para["N"]

  for i= 1:L

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
    for j =1:L

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


"""Generates the 1D hamiltonian MPO for the given system configuration"""
function init_ham(para::Dict, L::Int, disx::Vector{Float64}, disy::Vector{Float64}, sites; if_gate=false)
  # parameters 

  factor= para["TEBDfactor"]
  τ = para["τ"]
  QE = para["QE"]
  QN = para["QN"]
  headoverride = para["headoverride"]
  
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
  res = add_ee!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)
  res = add_ne!(res, para, L, disx, disy, sites, if_gate=if_gate, head=head, factor= factor, τ=τ)


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
    append!(res, reverse(res))
    return res

  end 
end

"""Generates a 2D system hamiltonian. Has the option to use the built-in lattice setup or manual setup"""
function init_ham(para::Dict, L::Vector{Int}, disx::Vector{Float64}, disy::Vector{Float64}, sites)
  # parameters 
  Lx = L[1]
  Ly = L[2]
  t = para["t"]
  λ_ee = para["int_ee"]
  λ_ne = para["int_ne"]
  ζ_ne, ζ_ee = para["ζ"]
  exch = para["exch"]
  decay = para["decay"]
  self = para["self_nuc"]
  manual = para["manual"]
  range = para["range"]
  QE = para["QE"]
  xscale = para["xscale"]
  CN = para["CN"]

  if QE > 0
    println("QE currently not supported on 2D sys")
  end 

  if Lx > Ly
    println("By set up Lx has to be smaller")
    exit()
  end 

  disx = reshape(disx, (Lx, Ly))
  disy = reshape(disy, (Lx, Ly))
  L = Lx * Ly

  λ_ne = λ_ne * CN / L
  ampo = OpSum()

  if !manual
    lattice = square_lattice(Lx, Ly; yperiodic = false)

    # important, lattice is a collection of bonds, so horizontal bonds: # = ( Nx - 1) * Ny, vertical bonds: # = (Ny - 1) * Nx
    # when we add terms to the DMRG H, we add them individually by the bond

    # hopping part
    # we use the bonds to encode hopping
    for b in lattice

      x1 = trunc(Int, b.x1)
      x2 = trunc(Int, b.x2)
      y1 = trunc(Int, b.y1)
      y2 = trunc(Int, b.y2)

      r = dis(x1, y1, x2, y2, disx, disy, xscale)
      hop = hopping(decay, r)

      #println(hop)
      # Hopping
      ampo += -t * hop, "C",b.s1,"Cdag",b.s2
      ampo += -t * hop, "C",b.s2,"Cdag",b.s1
    end

    # for the Coulombb interaction, it's more straightforward to use the number of the sites directly
    for j=1:L
      # E-E
      x1 = div(j - 1, Ly) + 1
      y1 = mod(j - 1, Ly) + 1

      for k=max(1, j- range):j-1

        x2 = div(k - 1, Ly) + 1
        y2 = mod(k - 1, Ly) + 1
        # delta function setting up the exchange between nearest neighbor
        ifexch = (1 -  ==(1, abs(x1 - x2) + abs(y1 - y2)) * exch )
        
        # add the ee interaction term one by one
        ampo += λ_ee * ifexch / ( dis(x1, y1, x2, y2, disx, disy, xscale) + ζ_ee),"N",j,"N",k

        
      end
      
      # N-E

      # check if self self_nuc

      if self
        cursum = 0
      else
        cursum = -λ_ne / ζ_ne
      end

      for l=1:L

        x2 = div(l - 1, Ly) + 1
        y2 = mod(l - 1, Ly) + 1
        cursum += λ_ne / ( dis(x1, y1, x2, y2, disx, disy, xscale) + ζ_ne)
      end

      ampo += -cursum, "N", j
    end

  else 

    # we always hop down along the long side so that to minimize the range of the hopping
    # by setup Lx <= Ly
    # hopping terms

    for j=1:L-1

      r = dis(j, j + 1, disx, disy)
      hop = hopping(decay, r)
  
      #println(hop)
      # Hopping
      ampo += -t * hop, "C",j,"Cdag",j+1
      ampo += -t * hop, "C",j+1,"Cdag",j
    end
  
    for j=1:L
      # E-E
      for k=1:j-1
  
        # delta function setting up the exchange between nearest neighbor
        ifexch = (1 -  ==(1, abs(j - k)) * exch )
        
        # add the ee interaction term one by one
        ampo += 2 * λ_ee * ifexch / ( dis(j, k, disx, disy) + ζ_ee),"N",j,"N",k
      end
      
      # N-E
  
      # check if self self_nuc
  
      if self
        cursum = 0
      else
        cursum = -λ_ne / ζ_ne
      end
  
      for l=1:L
        cursum += λ_ne / ( dis(j, l, disx, disy) + ζ_ne)
      end
  
      ampo += -cursum, "N", j
    end

  end 

  H = MPO(ampo,sites)
  return H
end
