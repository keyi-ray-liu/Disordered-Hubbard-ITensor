

using ITensors
using ITensorTDVP
using Random
using Test
using StatsBase
using ITensors.HDF5

# @testset "DMRG-X" begin

#   function HB(L)
#     ampo = OpSum()
#     hop = 1.0
#     t = 1.0

#     for j=1  :L-1
        
#         #println(hop)
#         # Hopping
#         ampo += -t * hop, "C",j,"Cdag",j+1
#         ampo += -t * hop, "C",j+1,"Cdag",j 
#     end

#     return ampo
#   end 

#   L = 12
#   N = 6

#   s = siteinds("Fermion", L, conserve_qns=true)
#   H = MPO(HB(L), s)

#   initstate = append!([ "Occ" for n=1:N] , ["Emp" for n=1:L -N])
#   ψ = randomMPS(s, initstate)

#   dmrg_x_kwargs = (
#     nsweeps=20, reverse_step=false, normalize=true, maxdim=20, cutoff=1e-10, outputlevel=0
#   )

#   ϕ = dmrg_x(ProjMPO(H), ψ; nsite=2, dmrg_x_kwargs...)

#   @test inner(ψ', H, ψ) / inner(ψ, ψ) ≈ inner(ϕ', H, ϕ) / inner(ϕ, ϕ) rtol = 1e-1
#   @test inner(H, ψ, H, ψ) ≉ inner(ψ', H, ψ)^2 rtol = 1e-7
#   @test inner(H, ϕ, H, ϕ) ≈ inner(ϕ', H, ϕ)^2 rtol = 1e-7

# end


# @testset "DMRG-X" begin
#   function heisenberg(n; h=zeros(n))
#     os = OpSum()
#     for j in 1:(n - 1)
#       os += 0.5, "S+", j, "S-", j + 1
#       os += 0.5, "S-", j, "S+", j + 1
#       os += "Sz", j, "Sz", j + 1
#     end
#     for j in 1:n
#       if h[j] ≠ 0
#         os -= h[j], "Sz", j
#       end
#     end
#     return os
#   end

#   n = 10
#   s = siteinds("S=1/2", n)

#   Random.seed!(12)

#   W = 12
#   # Random fields h ∈ [-W, W]
#   h = W * (2 * rand(n) .- 1)
#   H = MPO(heisenberg(n; h), s)

#   initstate = rand(["↑", "↓"], n)
#   ψ = MPS(s, initstate)

#   dmrg_x_kwargs = (
#     nsweeps=20, reverse_step=false, normalize=true, maxdim=20, cutoff=1e-10, outputlevel=0
#   )

#   ϕ = dmrg_x(ProjMPO(H), ψ; nsite=2, dmrg_x_kwargs...)

#   @test inner(ψ', H, ψ) / inner(ψ, ψ) ≈ inner(ϕ', H, ϕ) / inner(ϕ, ϕ) rtol = 1e-1
#   @test inner(H, ψ, H, ψ) ≉ inner(ψ', H, ψ)^2 rtol = 1e-7
#   @test inner(H, ϕ, H, ϕ) ≈ inner(ϕ', H, ϕ)^2 rtol = 1e-7

#   ϕ̃ = dmrg_x(ProjMPO(H), ϕ; nsite=1, dmrg_x_kwargs...)

#   @test inner(ψ', H, ψ) / inner(ψ, ψ) ≈ inner(ϕ̃', H, ϕ̃) / inner(ϕ̃, ϕ̃) rtol = 1e-1
#   @test inner(H, ϕ̃, H, ϕ̃) ≈ inner(ϕ̃', H, ϕ̃)^2 rtol = 1e-5
#   # Sometimes broken, sometimes not
#   # @test abs(loginner(ϕ̃, ϕ) / n) ≈ 0.0 atol = 1e-6
# end

# nothing



# @testset "Test" begin

#   function HB(L)
#     ampo = OpSum()
#     hop = 1.0
#     t = 10.0

#     for j=1  :L-1
        
#         #println(hop)
#         # Hopping
#         ampo += -t * hop, "C",j,"Cdag",j+1
#         ampo += -t * hop, "C",j+1,"Cdag",j 
#     end

#     return ampo
#   end 

#   function get_gate(s, L, τ)
#     # Make gates (1,2),(2,3),(3,4),...
#     gates = ITensor[]
#     hop = 1.0
#     t = 10.0

#     for j in 1:(L - 1)
#       s1 = s[j]
#       s2 = s[j + 1]
#       hj =
#         - t * hop * op("C", s1) * op("Cdag", s2) +
#         - t * hop * op("C", s2) * op("Cdag", s1)
#       Gj = exp(-im * τ / 2 * hj)
#       push!(gates, Gj)
#     end
#     # Include gates in reverse order too
#     # (N,N-1),(N-1,N-2),...
#     append!(gates, reverse(gates))

#     return gates
#   end 

#   L = 20
#   N = 10
#   τ = 0.1
#   start = 0.0
#   fin = 2.0
   
#   s = siteinds("Fermion", L, conserve_qns=true)
#   H = MPO(HB(L), s)
#   G = get_gate(s, L, τ)

#   prefix = pwd()
#   initstate = append!([ "Occ" for n=1:N] , ["Emp" for n=1:L -N])
#   ψ = randomMPS(s, initstate)
  
#   cutoff = 1e-9
  

#   ψkrylov = ψ
#   ψtrotter = ψ


#   for dt in start:τ:fin

#     println("time: $dt")
#     #ψ1 = tdvp(H, ψ, -1.0im * τ;  nsweeps=20, cutoff, nsite=2)
#     ψkrylov = tdvp(H,  -im * τ, ψkrylov;  cutoff, nsite=2, time_step= -im * τ/2, normalize=true)

#     ψtrotter = apply(G, ψtrotter; cutoff)
#     normalize!(ψtrotter)
    
#     println("Krylov occ: ", expect(ψkrylov, "N"))
#     println("Trotter occ: ", expect(ψtrotter, "N"))

#     wf = h5open(  prefix * "TDVP" * string(dt) * ".h5", "w")
#     write(wf, "psi", ψkrylov)
#     close(wf)

#     wf = h5open(  prefix * "Trotter" * string(dt) * ".h5", "w")
#     write(wf, "psi", ψtrotter)
#     close(wf)

#   end 

# end


function testRNG()

  NR = 16
  @show StatsBase.sample(1:NR, div(NR, 2), replace = false)

end 

function benchmark()

  paras = setpara(;L=12, N=6, ex=2)
  main(paras)

end 


function boson_fermion_test()
  function f(n, source, bulk)
    if n > length(source) && n <= length(source) + bulk
      return "Electron"

    else
      return "Boson"
    end
  end 

  source = [2, 2]
  drain = [1, 1]
  bulk = 2



  sites = siteinds( n -> f(n, source, bulk), length(source) + length(drain) + bulk)

  h = OpSum()
  h += 1.0, "Adag", length(source), "Cup", length(source) + 1
  h += 1.0, "Cdagup", length(source) + 1, "A", length(source)

  H = MPO(h, sites)

  
end 

function test()

  N = 5
  s = siteinds("Fermion", N; conserve_qns =true)
  vecs = ones( (N, 2))

  #@show vecs

  psi = MPS(s, 1)


  for j=1:N
    psi[j] = ITensor(vecs[j,:],  s[j])
  end
  orthogonalize!(psi,1)
  normalize!(psi)

  return psi
end 

# tcds = rand(10, 100) +  1im* rand(10, 100)

# num, L = size(tcds)

# λ_ee = 2.0
# ζ=0.5

# EE_matrix = reduce(hcat, [ [λ_ee / (abs(j - k) + ζ) for j in 1:L] for k in 1:L])
# EE_matrix -= Diagonal(EE_matrix)
# M = zeros(L, L)
# for i in 1:num

#   tcd = tcds[i, :]

#   r1 = dot(tcd, EE_matrix, tcd)

#   r3 = transpose(tcd) * EE_matrix * tcd

#   r2 = 0
#   for j in 1:L
#     for k in 1:L

#       if j != k
#         r2 += tcd[j] * tcd[k] * λ_ee / (abs(j - k) + ζ)
#         M[j, k] = λ_ee / (abs(j - k) + ζ)

#       end
#     end
#   end

#   println(r1,  "    ", r2, "   ", r3)

# end 

# unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

# s_len = 128
# d_len = 128
# t =  1.0
# source_offset =0.25
# drain_offset = -0.25

# source_energies = [ (2 * t * cos( k * pi / (s_len + 1) )+ source_offset, k, 1) for k in 1:s_len] 
# drain_energies = [ (2 * t * cos( k * pi / (d_len + 1) ) + drain_offset, k, 2) for k in 1:d_len] 

# energies, ks, LR = unzip(sort( vcat(source_energies, drain_energies) ))
# print(energies, ks, LR)

# LR = [1,1, 0, 0, 1, 0, 1]

# findall( x-> x==1, LR)

testRNG()