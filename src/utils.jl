function get_type_dict(type)

  op_str = Dict(
    1 => type == "Boson" ? "0" : "Emp",
    2 => type == "Boson" ? "1" : type == "Fermion" ? "Occ" : "Up",
    3 => type == "Boson" ? "2" : "Dn",
    4 => type == "Boson" ? "3" : "UpDn"
  )

  return op_str

end 


"""Calculates hopping
Return the strength of hooping, baseline is 1 without disorder"""
function disorder_hopping(decay::Float64, r::Float64)
  
  return exp( -(r - 1) * decay)
end

"""generates the locations of the NN down the flattened array direction"""
function get_nn(L::Vector{Int}, t; snake=false, geometry ="linear", spec_hop_struct=Dict{Int, Dict{Int, Float64}()}())
  nns = Dict{Int, Dict{Int, Float64}}()

  # get the size of the whole system
  total = prod(L)

  for site in 1:total

    nn = Dict{Int, Float64}()
    if geometry == "linear"
      # for each dimension, generates a nearest neighbor down the flattened array, if possible

      # we attempt to always make the smallest hopping possible, by sorting the L array
      # in the snake geometry, the numbering is always continuous, whether the next largest hopping happens is determined by its position
      if L != sort(L)
        error("for optimal performance, the input dimensions must be sorted. adjust QE position accordingly")
      end 

      curstride = 1
      next = curstride 

      for dim in L
        
        next *= dim

        if !snake
          if site + curstride <= total && site % next != 0
            #append!(nn, (site + curstride, t))
            nn[site + curstride] =  t
          end 

        else

          newstride = total

          if curstride == 1
            newstride = 1

          elseif site % curstride !=0
            newstride = 2 * curstride - (site % curstride -  1) - site % curstride
          end 

          if site + newstride <= total
            #append!(nn, (site + newstride, t))
            nn[site + newstride] = t 
          end 
        end 

        curstride = next
      end 


    # currently only available for 1D chain
    elseif geometry == "2ndNN"

      if site < total 
        #append!(nn, (site + 1, t))
        nn[site + 1] = t
      end

      if site < total - 1
        #append!(nn, (site + 2, t))
        nn[site + 2] = t
      end 


    elseif geometry == "loop-pbc"

      if site < total
        #append!(nn, (site + 1,t ))
        nn[site + 1] = t

      elseif site == total
        #append!(nn, (1, t))
        nn[1] = t

      end 

    # chain dim will e hardcoreded for now 2x2 squares
    elseif geometry == "chain"

      if (total - 1) % 3 != 0 
        error(ArgumentError("Not a valid chain length"))
      end 

      if site < total && (site - 1) % 3 == 0

        #append!(nn, (site + 1, t))
        #append!(nn, (site + 2, t))
        nn[site + 1] = t
        nn[site + 2] = t

      elseif site < total
        #append!(nn, (3 * div(site - 1, 3) + 4, t))
        nn[ 3 * div(site - 1, 3) + 4] = t
      end 



    elseif geometry == "tworing"

      #for two rings, 1, 2 are the shared sites
      # for the first ring, start 1, 2, 3, .... N//2 + 1
      # for the second ring, start 1, 2, N//2 + 2, ..., N
      half = div(total, 2)
      if site == 1
        nn[2] = t

      elseif site == 2
        nn[3] = t
        nn[ half + 2] = t

      elseif site == half + 1 || site == total
        nn[1] = t

      else
        nn[site + 1]= t
      end


    else
      error(ArgumentError("undefined geometry"))

    end 

    nns[ site ] = nn
  end 

  for spe in keys(spec_hop_struct)
    to_override = spec_hop_struct[ spe]

    for new in keys(to_override)
      nns[spe][new] = to_override[new]
    end 
  end   

  print(nns)

  return nns
end 

"""Calculates distance between site i and site j given the geometry of L"""
function dis(i::Int, j::Int, L::Vector{Int}, scales::Vector{Float64}, args...)

  # preprocess for better consistent // operation
  disorder = args
  i -= 1
  j -= 1

  curstride = 1
  nextstride = L[1]
  r2 = 0

  for d = 1: max(length(L), length(disorder))

    # if enough terms in L, get coordinates
    if d <= length(L)
      ci = div( i % nextstride, curstride)
      cj = div( j % nextstride, curstride)

      diffsite = ci - cj

      curstride = nextstride
      if d < length(L)
        nextstride *= L[d + 1]
      end 

    else
      diffsite = 0
    end 

    # if enough terms in disorder, get disorder
    if d <= length(disorder)
      diffdis = disorder[d][i + 1] - disorder[d][ j + 1]
    else
      diffdis = 0
    end 

    scale = d <= length(scales) ? scales[d]  : 1.0
    
    r2 += ((diffsite + diffdis) * scale )^2
  end 

  return sqrt(r2)

end

"""Calculates distance to QE at location"""
function dis(i::Int, QEloc::Vector{Float64}, L::Vector{Int}, scales::Vector{Float64}, args...)
  
    # preprocess for better consistent // operation
    disorder = args
    i -= 1
  
    curstride = 1
    nextstride = L[1]
    r2 = 0
  
    for d = 1: max(length(L), length(disorder))
  
      # if enough terms in L, get coordinates
      if d <= length(L)
        ci = div( i % nextstride, curstride)
        
        diffsite = ci - (d <= length(QEloc) ? QEloc[d] : 0)
  
        curstride = nextstride
        if d < length(L)
          nextstride *= L[d + 1]
        end 
  
      else
        diffsite = 0
      end 
  
      # if enough terms in disorder, get disorder
      if d <= length(disorder)
        diffdis = disorder[d][i + 1] 
      else
        diffdis = 0
      end 
  
      scale = d <= length(scales) ? scales[d]  : 1.0
      
      r2 += ((diffsite + diffdis) * scale )^2
    end 

    return sqrt(r2)
end 


""" check max bond dim of MPS """
function checkmaxbond(ψ, step, cnt)

  # s = length(ψ)
  # maxbond = 0

  # for i = 1:s
  #   maxbond = max( maxbond, maximum(size(ψ[i])))
  # end 

  # return maxbond
  start = 1 + step * (cnt  - 1)
  fin = min( length(ψ), step * cnt)

  return maximum([ size(tensor) for tensor in ψ[start:fin]])
end 

""" defines X operator"""
function ITensors.op!(Op::ITensor, ::OpName"x", ::SiteType"Fermion", s::Index)
  Op[s' => 2, s => 1] = 1.0
  return Op[s' => 1, s => 2] = 1.0
end

"""push gate"""
function gatefy!(gates, factor, g, τ)

  G = exp( - im * τ / factor * g)
  push!(gates, G)
end 

"""Parsing utility that grabs the max dim string"""
function stringmaxproc(output, tol)
  maxsweep = 0
  key = "maxlinkdim="
  enekey = "energy="
  output = split(output)

  for seg in output

    if length(seg) > length(key) && seg[begin:length(key) ] == key
      maxsweep = max( maxsweep, parse(Int64, seg[ length(key) + 1: end]))

    end
  end

  final_energy = []

  for seg in reverse(output)

    if length(final_energy) >= 2

      diff = abs( final_energy[1] - final_energy[2])
      if diff < tol
        println("final energy converged")

      else
        println("Warning: final energy not converged, difference: $diff")

      end

      break
    end

    if length(seg) > length(enekey) && seg[begin:length(enekey)] == enekey
      append!(final_energy, parse(Float64, seg[ length(enekey) + 1: end]))
    end

  end
  return maxsweep
end

"""Wrapper function for the evaluation of the std of Hamiltonian"""
function variance(H::MPO, psi::MPS)
  @suppress begin 
    var = inner(H, psi, H, psi) - inner(psi', H, psi) ^ 2 
    return var
  end
end 

linsolvemeasure!(o::AbstractObserver; kwargs...) = nothing
linsolvecheckdone!(o::AbstractObserver; kwargs...) = false

function linsolvemeasure!(obs::DMRGObserver; kwargs...)
  half_sweep = kwargs[:half_sweep]
  b = kwargs[:bond]
  psi = kwargs[:psi]
  truncerr = truncerror(kwargs[:spec])

  if half_sweep == 2
    N = length(psi)

    if b == (N - 1)
      for o in ops(obs)
        push!(measurements(obs)[o], zeros(N))
      end
      push!(truncerrors(obs), 0.0)
    end

    # when sweeping left the orthogonality center is located
    # at site n=b after the bond update.
    # We want to measure at n=b+1 because there the tensor has been
    # already fully updated (by the right and left pass of the sweep).
    wf = psi[b] * psi[b + 1]
    measurelocalops!(obs, wf, b + 1)

    if b == 1
      measurelocalops!(obs, wf, b)
    end
    truncerr > truncerrors(obs)[end] && (truncerrors(obs)[end] = truncerr)
  end
end


function count_ele(config)
  return count( i-> (i ==2), config) + count( i-> (i == 3), config) + 2 * count( i-> (i == 4), config)
end 

function Ukj(k, j, N)
  return sqrt( 2/ (N + 1)) * sin(  j * k * pi / (N + 1))

end 


"""generates the source and drain. Depending on the inelastic mix basis set up, we might add additional sites corresponding to adding Δk to k_i. Current operating mode: keyword for "inelastic_mode": "standard" """
function sd_gen(sd_hop; kwargs...)

  #which = get(kwargs, :which, "source")

  type = get(kwargs, :type, "Fermion")
  NR = get(sd_hop, "NR", 128)
  inelastic = get(sd_hop, "inelastic", false)

  occ1 = 1 + (type == "Fermion" ? 0 : 1)
  occ2 = occ1 + 1

  if !inelastic
    config = [ x in StatsBase.sample(1:NR, div(NR, 2), replace = false) ? occ1 : occ2 for x in 1:NR]

  else

    inelastic_para = get(sd_hop, "inelastic_para", Dict())
    kcnt = get(inelastic_para, "kcnt", 0)
    
    config = [ x in StatsBase.sample(1:NR, div(NR, 2), replace = false) ? occ1 : occ2 for x in 1:NR * (kcnt + 1)]
  end 



  return config
end 

""" 
function that sort energies of the reservoirs, either from scratch (ks = [] and LR = []), or with given order of ks. The latter is used when you want to set up the reservoir that's obtained with a different chemical potential that the current one.
"""
function get_mix_energy(para; ks = [], LR = [])

  unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

  sd_hop = para[:sd_hop]


  s_len = para[:s_len]
  d_len = para[:d_len]
  t = para[:t]

  source_offset = sd_hop["source_offset"]
  drain_offset = sd_hop["drain_offset"]

  if length(ks) != length(LR)
    error("input dimension does not match")
  end 

  # no ordering to follow, new calculations
  if length(ks) != s_len + d_len
    source_vals = [ (2 * t * cos( k * pi / (s_len + 1) )+ source_offset, k, 1) for k in 1:s_len] 
    drain_vals = [ (2 * t * cos( k * pi / (d_len + 1) ) + drain_offset, k, 2) for k in 1:d_len] 
    
    energies, ks, LR= unzip(sort( vcat(source_vals, drain_vals) ))

    # this is where the site starts in the whole lattice
    zeropoint = searchsortedfirst( energies, 0)

  # ks is supplied, keep ordering
  else
    
    offset(x) = x == 1 ? source_offset : drain_offset
    lens = [ s_len, d_len]

    energies = [ 2 * t * cos( ks[i] * pi / (lens[LR[i]] + 1) ) + offset(LR[i]) for i in 1: (s_len + d_len)]

    zeropoint = -1
  end 

  return energies, ks, LR, zeropoint
end 


function load_JSON(location)

  #temp_string = read(location, String)
  data = JSON3.read(location, Dict{String, Any} )

  print(typeof(data))
  return data
end 



"""function that calculates partial trace of MPO. The leftend indicates where the left partial trace ends ( traced up to, include), default 0 (no left partial trace); similar for right end ( include), default to 1e4 (no right partial trace)"""
function partial_tr( M::MPO; leftend::Int =0, rightend::Int=10^4, plev::Pair{Int,Int}=0 => 1, tags::Pair=ts"" => ts"")


  if leftend > rightend 
    error("leftend cannot be greater or equal to rightend")
  end 

  if leftend >= rightend - 1
    return tr(M; plev=plev, tags=tags)
  end 

  sys = deepcopy(M[ leftend + 1: rightend - 1])
  L = R = 1
  N = length(M)

  for left = 1:leftend
    L *= M[left]
    L = tr(L; plev=plev, tags=tags)
  end 

  for right =rightend:N
    R *= M[right]
    R = tr(R; plev=plev, tags=tags)
  end

  sys[1] *= L
  sys[end] *= R

  return MPO(sys)
end

function partial_contract(ψ::MPS, leftend::Int, rightend::Int)
  if leftend > rightend 
    error("leftend cannot be greater than rightend")
  end 

  if leftend >= rightend - 1
    return inner(ψ, ψ)
  end 

  sys = deepcopy(ψ[ leftend + 1: rightend - 1])

  ψdag = prime(linkinds, dag(ψ))
  L = R = ITensor(1.)

  for left = 1:leftend

    L *= ψ[left]
    L *= ψdag[left]
    #L *= prime(linkinds, dag(ψ[left]))

  end 


  for right =rightend:length(ψ)

    R *= ψ[right]
    R *= ψdag[right]
    
  end

  sys[1] *= L
  sys[end] *= R

  return sys
end 