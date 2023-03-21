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

#testlin()
load()
#check()
