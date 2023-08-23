"""Calculates hopping
Return the strength of hooping, baseline is 1 without disorder"""
function disorder_hopping(decay::Float64, r::Float64)
  
  return exp( -(r - 1)/decay)
end

"""generates the locations of the NN down the flattened array direction"""
function get_nn(site::Int, L::Vector{Int}; snake=false, geometry ="linear")
  nn = Int[]

  # get the size of the whole system
  total = prod(L)

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
          append!(nn, site + curstride)
        end 

      else

        newstride = total

        if curstride == 1
          newstride = 1

        elseif site % curstride !=0
          newstride = 2 * curstride - (site % curstride -  1) - site % curstride
        end 

        if site + newstride <= total
          append!(nn, site + newstride)
        end 
      end 

      curstride = next
    end 


  # currently only available for 1D chain
  elseif geometry == "2ndNN"

    if site < total 
      append!(nn, site + 1)
    end

    if site < total - 1
      append!(nn, site + 2)
    end 


  elseif geometry == "loop-pbc"

    if site < total
      append!(nn, site + 1)

    elseif site == total
      append!(nn, 1)

    end 



  # chain dim will e hardcoreded for now 2x2 squares
  elseif geometry == "chain"

    if (total - 1) % 3 != 0 
      error(ArgumentError("Not a valid chain length"))
    end 

    if site < total && (site - 1) % 3 == 0

      append!(nn, site + 1)
      append!(nn, site + 2)

    elseif site < total
      append!(nn, 3 * div(site - 1, 3) + 4)
    end 


  else
    error(ArgumentError("undefined geometry"))

  end 

  return nn
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
function checkmaxbond(ψ)

  s = length(ψ)
  maxbond = 0

  for i = 1:s
    maxbond = max( maxbond, maximum(size(ψ[i])))
  end 

  return maxbond
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

"""
Load QE file
"""
function load_qe()

  prefix = getworkdir()

  if !isfile(prefix * "QE.h5") 
    error(ArgumentError("Please provide the eigen basis function"))
  end 

  if !isfile(prefix * "QEex")
    error(ArgumentError("Please provide the eigen basis energy"))
  end 

  if !isfile(prefix * "QEvar")
    error(ArgumentError("Please provide the eigen basis variance"))
  end 

  staticenergy = vec(readdlm( prefix * "QEex"))

  staticwf= MPS[]
  staticwffile = h5open( prefix * "QE.h5")
  staticvar = vec(readdlm( prefix * "QEvar"))

  for key in sort(keys(staticwffile), by= x-> parse(Int, x[4:end]))
    append!(staticwf, [read(staticwffile, key, MPS )])
  end 

  if length(staticwf) != length(staticenergy) || length(staticwf) != length(staticvar)
    error(ArgumentError("QE wf and energy or var length mismatch!"))
  end 

  return staticwf, staticenergy, staticvar

end 

"""
Load eigenvectors and eigenvalues for the static hamiltonian eigenstates. \n
Return staticenergy, staticwf, overlaps
"""
function load_eigen(ψ)

  prefix = getworkdir()
  staticwf, staticenergy, _ = load_qe()

  overlaps = [ inner(ψ', staticwf[i]) for i in eachindex(staticwf)]
  println( "overlaps:", overlaps)
  writedlm( prefix * "overlaps", overlaps)
  
  println( "overlap sum:", sum( abs2.(overlaps)))
  return staticenergy, staticwf, overlaps

end 

function load_tcd()

  prefix = getworkdir()
  tcd_key = "TCD"
  tcd_dict = Dict()

  for file in filter(x->occursin(tcd_key,x), readdir(prefix))

    keys = split(file, "_")
    i = parse(Int, keys[end - 1])
    j = parse(Int, keys[end])
    
    tcd_dict[ (i, j) ] = readdlm( prefix * file)
  end 

  print(keys(tcd_dict))
  return tcd_dict

end 


