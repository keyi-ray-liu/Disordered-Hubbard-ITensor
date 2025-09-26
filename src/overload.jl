


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