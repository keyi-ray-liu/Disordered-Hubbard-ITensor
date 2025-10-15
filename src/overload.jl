


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








"""
    blocksites(psi, sites; range::UnitRange)

Fuse (block) the physical sites in `range` into a single supersite.
Returns (psiB, sitesB, map_unblock) where:
  - psiB is the new MPS,
  - sitesB are the new site indices (with one supersite),
  - map_unblock is an ITensor isometry you can use to expand the supersite back if needed.
"""
function blocksites(psi::MPS, sites; range::UnitRange)
  i, j = first(range), last(range)
  @assert 1 ≤ i ≤ j ≤ length(psi)

  # Move the orthogonality center to the block to keep gauge tidy
  orthogonalize!(psi, i)

  # Build a combiner that maps (s_i ⊗ s_{i+1} ⊗ … ⊗ s_j) → sB
  C = combiner((sites[k] for k in i:j)...; tags="Supersite($(i):$(j))")
  sB = commonind(C)  # the fused physical index

  # Contract the block of site tensors into one tensor T_block with left/right links and fused sB
  T_block = ITensor()
  for k in i:j
    T_block = (k == i ? psi[k] : T_block * psi[k])
  end
  # Replace the multiple site indices with the fused supersite using C'
  # (C' maps sB → s_i ⊗ … ⊗ s_j, so multiplying by C gives the fused index)
  T_block = T_block * C

  # Assemble the new MPS tensors:
  new_tensors = Any[]
  append!(new_tensors, psi[1:i-1])

  # Build the blocked site tensor with the correct left/right link structure
  # Identify the outer left/right links bounding the block
  left_link  = (i == 1      ? nothing : commonind(psi[i], psi[i-1]))  # L_{i-1}
  right_link = (j == length(psi) ? nothing : commonind(psi[j], psi[j+1]))  # L_j

  if left_link === nothing && right_link === nothing
    # Whole chain is one block: just a single-site MPS tensor
    push!(new_tensors, T_block)
  elseif left_link === nothing
    # Block at the beginning: order as (sB, right_link)
    T_block = permute(T_block, (sB, right_link))
    push!(new_tensors, T_block)
  elseif right_link === nothing
    # Block at the end: order as (left_link, sB)
    T_block = permute(T_block, (left_link, sB))
    push!(new_tensors, T_block)
  else
    # Middle block: order as (left_link, sB, right_link)
    T_block = permute(T_block, (left_link, sB, right_link))
    push!(new_tensors, T_block)
  end

  append!(new_tensors, psi[j+1:end])

  # Build new site index array with the supersite in place of i…j
  sitesB = Index[]
  append!(sitesB, sites[1:i-1])
  push!(sitesB, sB)
  append!(sitesB, sites[j+1:end])

  psiB = MPS(new_tensors)

  # map_unblock lets you expand operators/states back if needed:
  # it is the isometry C' : sB → (s_i ⊗ … ⊗ s_j)
  map_unblock = dag(C)

  return psiB, sitesB, map_unblock
end