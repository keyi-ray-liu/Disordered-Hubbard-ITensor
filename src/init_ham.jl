"""Generates the 1D hamiltonian MPO for the given system configuration"""
function init_ham(para::Dict, L::Int, disx::Vector{Float64}, disy::Vector{Float64}, sites)
  # parameters 
  t = para["t"]
  λ_ee = para["int_ee"]
  λ_ne = para["int_ne"]
  ζ_ne, ζ_ee = para["ζ"]
  exch = para["exch"]
  decay = para["decay"]
  self = para["self_nuc"]
  QE = para["QE"]
  QN = para["QN"]
  N = para["N"]
  CN = para["CN"]
  QEen = para["QEen"]
  dp = para["dp"]
  
  # set e-e interaction range
  range = para["range"]
  # if QE > 0, then at least left emitter
  head = (QE > 0)
  ampo = OpSum()

  λ_ne = λ_ne * CN / L

  # testblock
  ###### ###### ###### ###### ###### ######

  # if QE == 1
  #   head = 0
  #   L += 1
  #   disx = vcat([0.0], disx)
  #   disy = vcat([0.0], disy)
  # end 

  ######## ###### ###### ###### ###### ###### ######

  # adjust the site based on if there are left emitter
  for j=1  :L-1 

    r = dis(j, j + 1, disx, disy)
    hop = hopping(decay, r)

    #println(hop)
    # Hopping
    ampo += -t * hop, "C",j + head,"Cdag",j+1 + head
    ampo += -t * hop, "C",j+1 + head,"Cdag",j + head
  end

  for j=1 :L 
    # E-E
    for k= max(1, j - range):j-1

      # delta function setting up the exchange between nearest neighbor
      ifexch = (1 -  ==(1, abs(j - k)) * exch )
      
      # add the ee interaction term one by one
      ampo += λ_ee * ifexch / ( dis(j, k, disx, disy) + ζ_ee),"N",j + head,"N",k + head
    end
    
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

    ampo += -cursum, "N", j + head
  end

  
  # QE part
  # left QE
  if QE > 0

    # diagonal energy term

    ampo += QEen, "N", 1
    
    for left = 1 : L 
      
      r = dis(left, disx, disy)
      for all = 1 : L

        # r0 determines the overall 'weight' of sites
        r0 = dis(all, disx, disy)
        
        #off-diagonal two level transition term
        ampo  += - dp * r0 / r ^ 3, "x", 1, "N", left + 1, "N", all + 1
        
      end 
      
    end 

  end 

  # right QE
  if QE > 1
    ampo += QEen, "N", L + 2

    for right = 1 : L
      

      r = dis(right, disx, disy)

      for all = 1: L

        r0 = dis(right, disx, disy)
        ampo += - dp * (L + 1 - r0) / ( L + 1 - r) ^ 3, "x", L + 2, "N", right + 1, "N", all + 1
        #ampo += - dp / ( L + 1 - r) ^ 3, "Cdag", L + 2, "N", right + 1, "N", all + 1

      end 
    end 

  end 


  # penalty terms for non-QN conserving hamiltonian
  # in the form of (sum_i n_i - N) ^2
  if !QN
    Λ = 30.0

    for s1= 1:L
      # linear terms
      ampo += - 2 * Λ * N, "N", s1 + head

      # quadratic terms
      for s2 =1:L
        ampo += Λ, "N", s1 + head, "N", s2 + head
      end 

    end 

  end 


  H = MPO(ampo, sites)

  return H
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
