"""Calculates hopping
Return the strength of hooping, baseline is 1 without disorder"""
function hopping(decay::Float64, r::Float64)
  
  return exp( -(r - 1)/decay)
end

"""Calculates 1D distance"""
function dis(i::Int, j::Int, disx, disy)

  return sqrt( ( disx[i] + i - disx[j] - j) ^ 2 + (disy[i] - disy[j]) ^ 2 )

end

"""Calculates 2D distance"""
function dis(x1::Int, y1::Int, x2::Int, y2::Int, disx, disy)

  
  x = disx[x1, y1] + x1 - disx[x2, y2] - x2
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