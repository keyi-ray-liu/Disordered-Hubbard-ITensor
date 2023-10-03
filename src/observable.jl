"""Calculates correlation matrix between two wf"""
function cal_obs(ψ1, ψ2)

  # ψ1 is left, ψ2 is right

  L = length(ψ1)
  # type? fix if necessary
  corr = zeros( (L, L) )
  # get site indices

  
  s1 = siteinds(ψ1)
  s2 = siteinds(ψ2)

  for j = 1:L
    replaceind!(ψ1[j], s1[j], s2[j])
  end 

  # init post-operated state array
  allleft = []
  allright = []

  operator = "C"
  for site in 1:L

    ψleft = deepcopy(ψ1)
    ψright = deepcopy(ψ2)

    for prev in 1:site-1

      prevleft  = op("F", s2, prev) * ψleft[prev]
      ψleft[prev] = prevleft

      prevright = op("F", s2, prev) * ψright[prev]
      ψright[prev] = prevright

    end 
    
    
    curleft =  op( operator, s2, site) * ψleft[site]
    curright =  op( operator, s2, site) * ψright[site]
    #noprime!(curleft)
  
    ψleft[site] = curleft
    ψright[site] = curright

    #operator = "F * $operator "
    append!(allleft, [ψleft])
    append!(allright, [ψright])
    
  end 

  
  # explicitly calculate inner product of post-operated states

  for i in 1:L
    for j in 1:L
      corr[ i, j] = inner( allleft[i]', allright[j])
    end 
  end 


  return corr
end 

"""Calculates transition charge density between two wf"""
function cal_tcd(ψ1, ψ2; temp=true)

  # ψ1 is left, ψ2 is right

  L = length(ψ1)
  # type? fix if necessary
  tcd = zeros( ComplexF64, (L) )
  # get site indices

  #s1 = siteinds(ψ1)
  s2 = siteinds(ψ2)

  # for j = 1:L
  #   replaceind!(ψ1[j], s1[j], s2[j])
  # end 

  operator = "N"

  if !temp

    ψdag = dag(ψ1)
    
    
    lefts = []
    rights = []

    left = OneITensor()
    right = OneITensor()

    for i = 1:L

      append!(lefts, [left])
      append!(rights, [right])

      left = left *  ψdag[i] * ψ2[i]
      right = right * ψdag[ end + 1 - i] * ψ2[ end + 1 - i]

    end 

    reverse!(rights)
  end 


  # explicitly calculate inner product of post-operated states

  for i in 1:L
    #tcd[ i] = inner( ψ1', op( operator, s2, i), ψ2)

    Op = op( operator, s2, i)
    Nψ2 = apply( Op, ψ2)

    #tcd[i] = (lefts[i] * ψdag[i] * cur * rights[i])[]
    tcd[i] = inner( ψ1', Nψ2)

  end 

  #print(tcd)

  return tcd
end 


function time_corr_plot(paras)

  workdir = getworkdir()

  op1 = paras["op1"]
  op2 = paras["op2"]
  tag = paras["tag"]
  wftag = paras["wftag"]

  files = glob( wftag * "*.h5", workdir)

  if length(files) == 0
    error(ArgumentError("no time evolved wf found"))
  end 

  get_time(x) = parse(Float64, x[ length(workdir) + length(wftag ) + 1:end-3])
  sort!(files, by = x -> get_time(x))

  for file in files

    println("calculating file", file)
    wf = h5open(file, "r")
    ψ = read(wf, "psi", MPS)

    NN_corr = abs.(correlation_matrix(ψ, op1, op2))
    writedlm( workdir * tag * "_corr" * string(get_time(file)), NN_corr )
  end 


end 



function cal_current(ψ, para)

  N = para["N"]
  v = para["v"]

  Cupup = correlation_matrix(ψ, "Cdagup", "Cup")
  Cupdn = correlation_matrix(ψ, "Cdagup", "Cdn")
  Cdndn = correlation_matrix(ψ, "Cdagdn", "Cdn")

  exp_upup = imag(Cupup[1:N, N + 1])
  exp_dndn = imag(Cdndn[1:N, N + 1])
  exp_updn = imag(Cupdn[1:N, N + 1])
  exp_dnup = imag(conj.(transpose(Cupdn))[1:N, N + 1])

  U_matrix = reduce(hcat, [[ Ukj(k, j, N) for j in 1:N ] for k in 1:N])
  U_k1 = [ Ukj(k, 1, N) for k in 1:N]

  #println([ dot(U_k1, U_matrix, vec) for vec in (exp_upup, exp_dndn, exp_updn, exp_dnup)])
  I = - 2 * v * sum([ transpose(U_k1) * U_matrix * vec for vec in (exp_upup, exp_dndn, exp_updn, exp_dnup)])
  return I

end 

function cal_gpi(tcd_gs::Matrix, λ_ee, ζ, section)

  num, L = size(tcd_gs)
  gpi = zeros(ComplexF64, num, section)
  
  if typeof(section) == Int

    gap = div(L, section)
    section = [ ( 1 + 2 + gap * (i - 1), 2 + min(L - 4, gap * i)) for i in 1:section]

  end 
  
  print(section)

  EE_matrix = reduce(hcat, [ [λ_ee / (abs(j - k) + ζ) for j in 1:L] for k in 1:L])
  EE_matrix -= Diagonal(EE_matrix)


  for i in 1:num

    tcd = tcd_gs[i,:]
    gpi[i, :] = [transpose(tcd[start:fin]) * EE_matrix[start:fin, start:fin] * tcd[start:fin] for (start, fin) in section]
  end 

  return gpi
end 

function entropy_von_neumann(psi::MPS, b::Int)
  s = siteinds(psi)  
  orthogonalize!(psi, b)
  _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
  SvN = 0.0
  for n in 1:dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log(p)
  end
  return SvN
end

function cal_ee(ψ::MPS, system)

  L = length(ψ)

  soi = []

  if system == "QE"
    append!(soi, 2)
    append!(soi, 3)

    step = 10
    seg = div( L -4, step)

    for s in 1:seg 
      append!(soi, 2 + step * s)
    end 

    append!(soi, L - 1)

  elseif system == "SD"

    append!(soi, div(L, 2) + 1)

  end 

  return [ entropy_von_neumann(ψ, i ) for i in soi], [ maximum(size(ψ[ j])) for j in soi], soi

end 