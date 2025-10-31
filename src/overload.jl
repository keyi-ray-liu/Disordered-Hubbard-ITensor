


# function expect(psi::MPS, opname::String)
#   return expect(psi, opname)
# end 



# function myexpect(psi::MPS, selectors::Tuple{Dict}; sites=1:length(psi), site_range=nothing)
#   psi = copy(psi)
#   N = length(psi)
#   ElT = scalartype(psi)
#   s = siteinds(psi)

#   if !isnothing(site_range)
#     @warn "The `site_range` keyword arg. to `expect` is deprecated: use the keyword `sites` instead"
#     sites = site_range
#   end

#   site_range = (sites isa AbstractRange) ? sites : collect(sites)
#   Ns = length(site_range)
#   start_site = first(site_range)

#   #el_types = map(o -> ishermitian(op(o, s[start_site])) ? real(ElT) : ElT, selectors)

#   el_types = map(o -> ishermitian(op(o[1], s[start_site])) ? real(ElT) : ElT, selectors)

#   psi = orthogonalize(psi, start_site)
#   norm2_psi = norm(psi)^2
#   iszero(norm2_psi) && error("MPS has zero norm in function `expect`")

#   ex = map((o, el_t) -> zeros(el_t, Ns), selectors, el_types)
#   for (entry, j) in enumerate(site_range)
#     psi = orthogonalize(psi, j)
#     for (n, selector) in enumerate(selectors)

#         opname = selector[j]

#         #oⱼ = adapt(datatype(psi[j]), op(opname, s[j]))
#         #val = inner(psi[j], apply(oⱼ, psi[j])) / norm2_psi
#         o = op(opname, s[j])
#         val = inner(psi[j], apply(o, psi[j])) / norm2_psi
#         ex[n][entry] = (el_types[n] <: Real) ? real(val) : val
#     end
#   end

#   if sites isa Number
#     return map(arr -> arr[1], ex)
#   end
#   return ex
# end
function myexpect(psi::MPS, op::String, replacement::String; kwargs...)
  return first(myexpect(psi, (op,), replacement; kwargs...))
end


function myexpect(psi::MPS, ops, replacement; sites=1:length(psi), site_range=nothing)
  psi = copy(psi)
  N = length(psi)
  ElT = scalartype(psi)
  s = siteinds(psi)

  if !isnothing(site_range)
    @warn "The `site_range` keyword arg. to `expect` is deprecated: use the keyword `sites` instead"
    sites = site_range
  end

  site_range = (sites isa AbstractRange) ? sites : collect(sites)
  Ns = length(site_range)
  start_site = first(site_range)

  el_types = map(o -> ishermitian(op(o, s[start_site])) ? real(ElT) : ElT, ops)

  psi = orthogonalize(psi, start_site)
  norm2_psi = norm(psi)^2
  iszero(norm2_psi) && error("MPS has zero norm in function `expect`")

  ex = map((o, el_t) -> zeros(el_t, Ns), ops, el_types)
  for (entry, j) in enumerate(site_range)
    psi = orthogonalize(psi, j)
    for (n, opname) in enumerate(ops)

      val = nothing 

      try
        oⱼ = op(opname, s[j])
        val = inner(psi[j], apply(oⱼ, psi[j])) / norm2_psi

      catch e
        oⱼ = op(replacement, s[j])
        val = inner(psi[j], apply(oⱼ, psi[j])) / norm2_psi
      end 

      ex[n][entry] = (el_types[n] <: Real) ? real(val) : val
    end
  end

  if sites isa Number
    return map(arr -> arr[1], ex)
  end
  return ex
end



function ITensors.space(::SiteType"TLS";
                        conserve_qns=false)
  if conserve_qns
    return [QN("N",0)=>2]
  end
  return 2
end

ITensors.op(::OpName"N",::SiteType"TLS") =
  [1 0; 0 0]

ITensors.op(::OpName"T",::SiteType"TLS") =
  [0 1; 1 0]


ITensors.state(::StateName"Emp", ::SiteType"TLS") = [1.0, 0.0]








# Given: psi (MPS), sites (site indices), and a vector 'target_order' of site labels
# Step 1: build the list of adjacent swaps to turn current order into target_order
function adjacent_swaps_from_permutation(curr::Vector{Int}, target::Vector{Int})
    swaps = Tuple{Int,Int}[]
    pos = Dict(v => i for (i,v) in enumerate(curr))
    for t in target
        i = pos[t]
        while i > 1 && curr[i-1] != target[i-1]
            # swap positions (i-1,i)
            push!(swaps, (i-1,i))
            curr[i-1], curr[i] = curr[i], curr[i-1]
            pos[curr[i-1]] = i-1; pos[curr[i]] = i
            i -= 1
        end
    end
    return swaps
end

# Example SWAP gate for spin-1/2 using (I+σ·σ)/2
function spin_swap_gate(sites, k)
    s1 = sites[k]; s2 = sites[k+1]
    I1 = op("Id", s1);  I2 = op("Id", s2)
    X1 = op("Sx", s1);  X2 = op("Sx", s2)
    Y1 = op("Sy", s1);  Y2 = op("Sy", s2)
    Z1 = op("Sz", s1);  Z2 = op("Sz", s2)
    G  = 0.5 * (I1 * I2 + X1 * X2 + Y1 * Y2 + Z1 * Z2)  # two-site ITensor
    return G
end

# ---------- Spinless fermion FSWAP on sites k,k+1 ----------
# Basis order for each site: |0>, |1>
function fswap_gate_fermion(sites, k)
    s1 = sites[k]
    s2 = sites[k+1]
    s1p = prime(s1)
    s2p = prime(s2)
    G = ITensor(s1p, s2p, dag(s1), dag(s2))

    # Helper to get state index
    i0_1 = stateind(s1, "0"); i1_1 = stateind(s1, "1")
    i0_2 = stateind(s2, "0"); i1_2 = stateind(s2, "1")
    o0_1 = stateind(s1p, "0"); o1_1 = stateind(s1p, "1")
    o0_2 = stateind(s2p, "0"); o1_2 = stateind(s2p, "1")

    # |00> -> |00|
    G[o0_1, o0_2, i0_1, i0_2] = 1.0
    # |01> -> |10|
    G[o1_1, o0_2, i0_1, i1_2] = 1.0
    # |10> -> |01|
    G[o0_1, o1_2, i1_1, i0_2] = 1.0
    # |11> -> -|11|
    G[o1_1, o1_2, i1_1, i1_2] = -1.0

    return G
end

# ---------- Spinful electron FSWAP on sites k,k+1 ----------
# States: "Emp"(even), "Up"(odd), "Dn"(odd), "UpDn"(even)
# Sign = (-1)^(parity(a)*parity(b)) where parity=1 for odd, 0 for even.
function fswap_gate_electron(sites, k)
    s1 = sites[k]
    s2 = sites[k+1]
    s1p = prime(s1)
    s2p = prime(s2)
    G = ITensor(s1p, s2p, dag(s1), dag(s2))

    labels = ("Emp","Up","Dn","UpDn")
    parity = Dict("Emp"=>0, "Up"=>1, "Dn"=>1, "UpDn"=>0)

    for a in labels, b in labels
        # input |a⟩_k |b⟩_{k+1}  -> output |b⟩_k |a⟩_{k+1} * sign
        sign = ((parity[a]*parity[b]) % 2 == 1) ? -1.0 : 1.0
        ia1 = stateind(s1, a); ib2 = stateind(s2, b)
        ob1 = stateind(s1p, b); oa2 = stateind(s2p, a)
        G[ob1, oa2, ia1, ib2] = sign
    end
    return G
end

# Fermionic FSWAP: build explicitly in the local basis or using parity-aware ops.
# (Sketch — details depend on your site type; core idea is to fill the 4-index tensor.)

function apply_swap_chain!(psi::MPS,  swaps::Vector{Tuple{Int,Int}})

    sites = siteinds(psi)

    for (k, k1) in swaps
        @assert k1 == k+1

        
        # Move orthogonality center near k
        psi = orthogonalize(psi, k)

        # Build the proper two-site gate:
        G = fswap_gate_fermion(sites, k) 

        # Apply gate EXACTLY (no truncation)
        psi = apply(G, psi; cutoff=0.0)  # allow temporary bond growth

        # Optionally: re-orthogonalize here or at the end
        psi = orthogonalize(psi, k)
    end
    return psi
end


function swapstate(curr, target, psi)

  swaps = adjacent_swaps_from_permutation(copy(curr), target)
  psi_perm = apply_swap_chain!(psi,  swaps)

  # Single final compression if bonds exploded:
  psi_perm = compress(psi_perm; cutoff=1e-12, maxdim=your_cap)

  return psi_perm
end 