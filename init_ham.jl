"""return the sites indices for further use"""
function init_site(para::Dict)
  L = para["L"]

  if typeof(L) == Int
    sites = siteinds("Fermion", L, conserve_qns=true)

  else
    sites = siteinds("Fermion", L[1] * L[2], conserve_qns=true)

  end 
  return sites
end 



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

  # set e-e interaction range
  range = para["range"]

  ampo = OpSum()
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
    for k= max(1, j - range):j-1

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

  H = MPO(ampo,sites)

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

  if Lx > Ly
    println("By set up Lx has to be smaller")
    exit()
  end 

  disx = reshape(disx, (Lx, Ly))
  disy = reshape(disy, (Lx, Ly))
  L = Lx * Ly

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

      r = dis(x1, y1, x2, y2, disx, disy)
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
        ampo += 2 * λ_ee * ifexch / ( dis(x1, y1, x2, y2, disx, disy) + ζ_ee),"N",j,"N",k

        
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
        cursum += λ_ne / ( dis(x1, y1, x2, y2, disx, disy) + ζ_ne)
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
