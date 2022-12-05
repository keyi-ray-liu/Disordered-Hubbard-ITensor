"""Calculates hopping
Return the strength of hooping, baseline is 1 without disorder"""
function hopping(decay::Float64, r::Float64)
  
  return exp( -(r - 1)/decay)
end

"""Calculates 1D distance"""
function dis(i::Int, j::Int, disx, disy)

  return sqrt( ( disx[i] + i - disx[j] - j) ^ 2 + (disy[i] - disy[j]) ^ 2 )

end

"""Calculates 1D distance to left QE"""
function dis(i::Int, disx, disy)
  return sqrt( ( disx[i] + i + 1) ^ 2 + (disy[i]) ^ 2)
end 

"""Calculates 2D distance, with potential scaling on x"""
function dis(x1::Int, y1::Int, x2::Int, y2::Int, disx, disy, xscale)

  
  x = disx[x1, y1]- disx[x2, y2]  + (x1  - x2) * xscale
  y = disy[x1, y1] + y1 - disy[x2, y2] - y2

  
  return sqrt( x^2 + y^2 )

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
  var = inner(H, psi, H, psi) - inner(psi, H, psi) ^ 2
  return sqrt( abs(var))

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