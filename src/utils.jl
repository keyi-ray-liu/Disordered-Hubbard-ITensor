function get_type_dict(type)

    op_str = Dict(
      1 => type == "Boson" ? "0" : "Emp",
      2 => type == "Boson" ? "1" : type == "Fermion" ? "Occ" : "Up",
      3 => type == "Boson" ? "2" : "Dn",
      4 => type == "Boson" ? "3" : "UpDn"
    )
  
    return op_str
  
end 

"""Set the work directory"""
function getworkdir()

  workdir = pwd() * "/work/"

  if !isdir(workdir)
    mkdir(workdir)
  end 

  return workdir

end 


load_JSON(location) = JSON3.read(location, Dict{String, Any} )

"""Wrapper function for the evaluation of the std of Hamiltonian"""
function variance(H::MPO, psi::MPS)
  @suppress begin 
    var = inner(H, psi, H, psi) - inner(psi', H, psi) ^ 2 
    return var
  end
end


function load_ψ(output::String; tag ="psi1")
  workdir = getworkdir()

  if length(output) > 3 && output[end-2:end] == ".h5"
    wf = h5open(workdir * output, "r")
  else
    wf = h5open( workdir * output * ".h5", "r")
  end 

  ψ = read(wf, tag, MPS)

  return ψ
end 

function load_ψ(t::Float64; tag ="psi1")
  workdir = getworkdir()
  wf = h5open( workdir * "tTDVP" * string(t) * ".h5", "r")
  ψ = read(wf, tag, MPS)

  return ψ
end 

function load_plsmon(output)
  ex = readdlm( getworkdir() * output * "ex")

  return ex[2] - ex[1]

end 

checkexist(output) = isfile( getworkdir() * output * ".h5")


# returns the precise site number of the ex QE site, beginning of the subsystem and end of subsystem, corresponding to site j
function get_sys_loc(sys::QE_flat_SIAM, j::Int) 

  if j <= left(sys)

    mod = (j - 1) % (siteseach(sys) + QESITES)
    qe_loc = j - mod 
    chain_begin = qe_loc + QESITES
    chain_end = chain_begin + siteseach(sys) - 1

  elseif j > left(sys) + 1

    mod = (j - 2) % (siteseach(sys) + QESITES)
    chain_begin = j - mod 
    chain_end = chain_begin + siteseach(sys) - 1
    qe_loc = chain_end + 1


  else
    qe_loc = chain_begin = chain_end = -1
  end 

  return qe_loc, chain_begin, chain_end

end 

function gen_graph(sys::QE_G_SIAM)

  g = NamedGraph()

  for i in 1:legleft(sys) + legright(sys)
    for j in 1:siteseach(sys) + QESITES
      add_vertex!(g, (i, j))
    end 
  end 

  center = (legleft(sys) + legright(sys) + 1, 1)
  add_vertex!(g, center)

  for i in 1:legright(sys) + legright(sys)

    add_edge!(g, (center, (i, 1)))
    for j in 1: siteseach(sys) +QESITES - 1
      add_edge!(g, ((i, j), (i, j + 1)))
    end 
  end 

  @show vertices(g)

  return g
  #@visualize g

end 

sitemap(sys::systems, j) = j
sitemap(sys::QE_G_SIAM, j) = sitemap(sys)[j]

"""maps the flattened index to graph index"""


function get_sitemap(sys::QE_flat_SIAM)

  d = Dict()

  for j in 1:get_systotal(sys)

    full = QESITES + siteseach(sys)

    if j < left(sys) + 1
  
      leg = div(j - 1, full) + 1
      idx = full - (j - 1) % full
  
    elseif j == left(sys) + 1
  
      leg = legleft(sys) + legright(sys) + 1
      idx = 1
  
    else
  
      j -= 1
      leg = div(j - 1, full) + 1
      idx = (j - 1) % full + 1
      j += 1
  
    end 
  
    d[(leg, idx)] = j
    d[j] = (leg, idx)


  end 

  @show d
  return d
end 


