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

    elseif geometry == "twosquare"

      each = div(total , 2)
      
      if site % each == 0
        #append!(nn, (1, t))
        nn[site - each + 1] = t

      else
        nn[site + 1] = t

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


function get_mix_energy(para)

  unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))
  sd_hop = para[:sd_hop]
  s_len = length(para[:source_config])
  d_len = length(para[:drain_config])
  t = para[:t]

  source_offset = sd_hop["source_offset"]
  drain_offset = sd_hop["drain_offset"]

  source_energies = [ (2 * t * cos( k * pi / (s_len + 1) )+ source_offset, k) for k in 1:s_len] 
  drain_energies = [ (2 * t * cos( k * pi / (d_len + 1) ) + drain_offset, k) for k in 1:d_len] 
  
  energies, ks = unzip(sort( vcat(source_energies, drain_energies) ))
  return energies, ks
end 