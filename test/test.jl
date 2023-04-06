using ITensors
using ITensorTDVP
using Test

function testlin()
  @testset "Linsolve" begin
    cutoff = 1E-11
    maxdim = 8
    nsweeps = 2

    N = 8
    s = siteinds("S=1/2", N)

    os = OpSum()
    for j in 1:(N - 1)
      os += 0.5, "S+", j, "S-", j + 1
      os += 0.5, "S-", j, "S+", j + 1
      os += "Sz", j, "Sz", j + 1
    end
    H = MPO(os, s)

    # Correct x is x_c
    x_c = randomMPS(s; linkdims=4)
    # Compute b
    b = apply(H, x_c; cutoff)

    x0 = randomMPS(s; linkdims=10)
    x = linsolve(H, b, x0; cutoff, maxdim, nsweeps, ishermitian=true, solver_tol=1E-6)

    #@show linkdims(x)
    @show norm(x - x_c)
    @test norm(x - x_c) < 1E-4
  end
end 

function load()

  baseline = h5open("QE.h5", "r")
  states = sort(keys(baseline), by= x -> parse(Int, x[4:end]))

  print(states)
  ex = length(states)
  exwf = Vector{MPS}(undef, ex)

  for (i, key) in enumerate(states)
    exwf[i] = read(baseline, key, MPS)
  end 

  s = length(exwf[15])

  for i=1:s
    println(size(exwf[15][s]))

  end 
end 

function dis(i::Int, QEoffset::Float64)
  return sqrt( (  i + QEoffset ) ^ 2 )
end 

function check()

  QEoffset = 1.0 
  QE = 2

  L = 50
  lefts = rights = []
  for left = 1:L
    r = dis(left, QEoffset)
    
    append!(lefts, r)
  end 

  for right = 1:L
    r  = ( L - 1 + QE * (QEoffset + 1) -  dis(right, QEoffset))  
    append!(rights, r)
  end 

  print( lefts == reverse!(rights))
end 

function get_nn(site::Int, L::Vector{Int})
  nn = Int[]

  # get the size of the whole system
  total = prod(L)
  # for each dimension, generates a nearest neighbor down the flattened array, if possible

  # we attempt to always make the smallest hopping possible, by sorting the L array
  Lsort = sort(L)

  curstride = 1
  next = curstride 

  for dim in Lsort
    
    next *= dim
    if site + curstride <= total && site % next != 0
      append!(nn, site + curstride)
    end 

    curstride = next
  end 

  println("$site, $nn")
  return nn
end 

function dis(i::Int, j::Int, L::Vector{Int}, scales, args...)

  # preprocess for better consistent operation
  disorder = args
  i -= 1
  j -= 1

 
  curstride = 1
  nextstride = L[1]
  r2 = 0

  for d = 1: max(length(L), length(disorder))

    # if enough terms in L, get coordinates
    if d <= length(L)
      ci = div( i % nextstride, curstride)
      cj = div( j % nextstride, curstride)

      diffsite = ci - cj

      curstride = nextstride
      if d < length(L)
        nextstride *= L[d + 1]
      end 
      
    else
      diffsite = 0
    end 

    # if enough terms in disorder, get disorder
    if d <= length(disorder)
      diffdis = disorder[d][i + 1 ] - disorder[d][ j + 1 ]
    else
      diffdis = 0
    end 

    scale = d <= length(scales) ? scales[d]  : 1.0
    
    r2 += ((diffsite + diffdis) * scale )^2
  end 


  return sqrt(r2)

end

function print_nn()

  L = [2, 30]
  t = prod(L)

  for i = 1: t
    get_nn(i, L)
  end 
end 
#testlin()
#load()
#check()

function print_dis()

  L = [2, 12]
  tot = prod(L)
  d = zeros(tot)
  for i = 1: tot
    for j =i + 1:  tot

      res = dis(i, j, L, [2.0], d, d)
      println("$i, $j, $res")

    end 
  end 


end 

function checknn()
  L = [15]
  t = prod(L)

  for i = 1:t
    for j =1:i-1

      println( "$i, $j, ", 1 - 0.2 * ( i in get_nn(j, L)))
    end 
  end 
end 

print_nn()
#print_dis()

#checknn()