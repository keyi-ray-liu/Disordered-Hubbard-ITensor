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