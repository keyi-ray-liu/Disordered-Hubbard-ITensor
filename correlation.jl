"""Calculates correlation matrix between two wf"""
function cal_corr(psi1, psi2)

  # psi1 is left, psi2 is right

  L = length(psi1)
  # type? fix if necessary
  corr = zeros( (L, L) )
  # get site indices

  s1 = siteinds(psi1)
  s2 = siteinds(psi2)

  for j = 1:L
    replaceind!(psi1[j], s1[j], s2[j])
  end 

  # init post-operated state array
  allleft = []
  allright = []

  operator = "C"
  for site in 1:L

    psileft = deepcopy(psi1)
    psiright = deepcopy(psi2)

    for prev in 1:site-1

      prevleft  = op("F", s2, prev) * psileft[prev]
      psileft[prev] = prevleft

      prevright = op("F", s2, prev) * psiright[prev]
      psiright[prev] = prevright

    end 
    
    
    curleft =  op( operator, s2, site) * psileft[site]
    curright =  op( operator, s2, site) * psiright[site]
    #noprime!(curleft)
  
    psileft[site] = curleft
    psiright[site] = curright

    #operator = "F * $operator "
    append!(allleft, [psileft])
    append!(allright, [psiright])
    
  end 

  
  # explicitly calculate inner product of post-operated states

  for i in 1:L
    for j in 1:L
      corr[ i, j] = inner( allleft[i]', allright[j])
    end 
  end 


  return corr
end 
