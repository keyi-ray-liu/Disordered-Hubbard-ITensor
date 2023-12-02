
"""add QN terms if necessary"""
function add_qn!(res, para::Dict,  L::Vector{Int}, sites; if_gate=false, head=0, factor=2, τ=0.1, Λ=30)
  
  println("Adding QN conserving")
  N = para["N"]
  Ltotal = prod(L)

  for i= 1:Ltotal

    p1 = i + head
    s1 = sites[p1]
    # linear terms

    if !if_gate
      res += - 2 * Λ * N, "N", s1 + head
    else
      li = - 2 * Λ * N * op( "N", s1)
      gatefy!(res, factor, li, τ)
    end 
      

    # quadratic terms
    for j =1:Ltotal

      p2 = j + head
      s2 = sites[p2]

      if !if_gate
        res += Λ, "N", p1, "N", p2
      else

        if s1 != s2
          lij = Λ * op("N", s1) * op("N", s2)
        end 
        gatefy!(res, factor, lij, τ)

      end 
    end 

  end 

  return res

end 

