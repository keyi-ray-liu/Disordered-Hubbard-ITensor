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