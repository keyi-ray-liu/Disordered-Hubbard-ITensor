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
  tcd = zeros( (L) )
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
    Op /= norm(Op)
    
    cur =  Op * ψ2[i] 
    noprime!(cur)

    new = copy(ψ2)
    new[i] = cur

    #tcd[i] = (lefts[i] * ψdag[i] * cur * rights[i])[]
    tcd[i] = inner( ψ1', new)

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

