
"""
position!(P::AbstractProjMPO, psi::MPS, phi::MPS, pos::Int)

This function deals with the situation where the left and right (bra and ket) of the MPO operator is not the same
left state is provided as psi, right state is provided as phi

Important distinction, we should return the result as an MPS instead of an MPO (that is, to be used as b in Ax = b, instead of A)

That is, instead of 
````
o--o--o-      -o--o--o--o--o--o <psi|
|  |  |  |  |  |  |  |  |  |  |
o--o--o--o--o--o--o--o--o--o--o H
|  |  |  |  |  |  |  |  |  |  |
o--o--o-      -o--o--o--o--o--o |psi>
````

We want 
````
o--o--o-      -o--o--o--o--o--o <psi|
|  |  |  |  |  |  |  |  |  |  |
o--o--o--o--o--o--o--o--o--o--o H
|  |  |  |  |  |  |  |  |  |  |
o--o--o--o--o--o--o--o--o--o--o |phi>
````
"""
nsite(P::AbstractProjMPO) = P.nsite

function position!(P::AbstractProjMPO, psi::MPS, phi::MPS, pos::Int)
  L = makeL!(P, psi, phi, pos - 1)
  R = makeR!(P, psi, phi, pos + nsite(P))

  # start of the left end
  if L == nothing
    return phi[pos] * phi[pos + 1] * R

  # start of the right end
  elseif R == nothing
    return L * phi[pos] * phi[pos + 1]
  
  # in the middle
  else
    return L * phi[pos] * phi[pos + 1] * R

  end 
end


function _makeL!(P::AbstractProjMPO, psi::MPS, phi::MPS, k::Int)
  # Save the last `L` that is made to help with caching
  # for DiskProjMPO
  ll = P.lpos

  # here we change the greater equal sign to greater, as additional operation has to be done on the k site to make it the correct dimensions
  if ll > k
    # Special case when nothing has to be done.
    # Still need to change the position if lproj is
    # being moved backward.
    P.lpos = k
    return nothing
  end
  # Make sure ll is at least 0 for the generic logic below
  ll = max(ll, 0)
  L = lproj(P)
  while ll < k
    L = L * phi[ll + 1] * P.H[ll + 1] * dag(prime(psi[ll + 1]))
    P.LR[ll + 1] = L
    ll += 1
  end

  # Needed when moving lproj backward.
  P.lpos = k
  return L
end

function makeL!(P::AbstractProjMPO, psi::MPS, phi::MPS, k::Int)
  L = _makeL!(P, psi, phi, k)
  #return P
  return L
end

function _makeR!(P::AbstractProjMPO, psi::MPS, phi::MPS, k::Int)
  # Save the last `R` that is made to help with caching
  # for DiskProjMPO
  rl = P.rpos
  if rl ≤ k
    # Special case when nothing has to be done.
    # Still need to change the position if rproj is
    # being moved backward.
    P.rpos = k
    return nothing
  end
  N = length(P.H)
  # Make sure rl is no bigger than `N + 1` for the generic logic below
  rl = min(rl, N + 1)
  R = rproj(P)
  while rl > k
    R = R * phi[rl - 1] * P.H[rl - 1] * dag(prime(psi[rl - 1]))
    P.LR[rl - 1] = R
    rl -= 1
  end
  P.rpos = k
  return R
end

function makeR!(P::AbstractProjMPO, psi::MPS, phi::MPS, k::Int)
  R = _makeR!(P, psi, phi, k)
  #return P
  return R
end
