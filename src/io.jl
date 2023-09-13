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

