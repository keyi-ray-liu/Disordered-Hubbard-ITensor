using Lazy
using DelimitedFiles
using Glob
using ITensors, ITensorMPS 
using JSON3
using LinearAlgebra
using StatsBase
using Suppressor
using ITensorGaussianMPS
using Random
using HDF5
using Logging
using Observers: observer
using StableRNGs: StableRNG
using Test: @test, @test_broken, @testset
using ITensorUnicodePlots: @visualize

include("systems.jl")
include("QEsystems.jl")
include("initial.jl")
include("Onsite.jl")
include("Hopping.jl")
include("Densitydensity.jl")
include("solve.jl")
include("simulation.jl")
include("Hubbard.jl")
include("QEterms.jl")
include("utils.jl")
include("observables.jl")
include("DPT.jl")
include("SD.jl")
include("QErun.jl")
include("NF.jl")
include("Chain.jl")
include("specific.jl")
include("basis.jl")
include("QEutil.jl")
include("test.jl")
include("SqChain.jl")

cd(@__DIR__)

using ITensors

function test_corr_MPO()
 

    # N = 4
    # s = siteinds("Electron", N; conserve_qns = true )
    # state = [ isodd(i) ? "Emp" : "Up" for i in 1:N]
    # psi1 = randomMPS(s, state)
 
    # state2 = [ iseven(i) ? "Emp" : "Up" for i in 1:N]
    # psi2 = randomMPS(s, state2)
 
    # psi = add( sqrt(0.7) * psi1, sqrt(0.3) * psi2)
    psi = load_ψ("wf.h5")
    s = siteinds(psi)
    N = length(s)
    @show expect(psi, "Sz")
 
    a = OpSum()
    for i in 1:N
        for j in 1:N
            a += 1.0, "Sz", i, "Sz", j
        end
    end
 
    A = MPO(a, s)
 
    @show inner(psi', A, psi)
    @show sum(correlation_matrix(psi, "Sz", "Sz"))
 
end



function total_spin_squared_mpo(sites)
    Num = length(sites)
    os = OpSum()

    for i in 1:Num
        # Single-site spin squared: Sx^2 + Sy^2 + Sz^2 = 3/4 when site is singly occupied
        # os += 0.75, "Id", i
        # os += 0.75, "Nup", i
        # os += 0.75, "Ndn", i
        # os += -1.5, "Nup", i, "Ndn", i


        # for j in i+1:N
        for j in 1:Num
            # os += 0.25, "Nup", i, "Nup", j
            # os += -0.25, "Ndn", i, "Nup", j
            # os += -0.25, "Nup", i, "Ndn", j
            # os += 0.25, "Ndn", i, "Ndn", j
            # # S_i · S_j = Sz_i Sz_j + 0.5*(S+_i S-_j + S-_i S+_j)
            os += 1.0, "Sz", i, "Sz", j
            os += 0.5, "S+", i, "S-", j
            os += 0.5, "S-", i, "S+", j
            # os += 1.0, "Sz", j, "Sz", i
            # # os += 0.5, "S+", j, "S-", i
            # # os += 0.5, "S-", j, "S+", i
            
        end
    end
    return MPO(os, sites)
end


@time begin

# test_corr_MPO()
let
    # Num = 16
    # sites = siteinds("Electron", Num)
    # psi = randomMPS(sites)  # Or your DMRG result
    psi = load_ψ("wf.h5")
    # psi = read("/Users/ynl42/Documents/GitHub/Disordered-Hubbard-ITensor/src/work/wf.h5", MPS)

    # @show expect(psi, "Nup")
    # @show expect(psi, "Sz")
    sites = siteinds(psi) 
    # @show length(sites)

    # corr = 0
    # @show corr += sum((correlation_matrix(psi, "Sz", "Sz")))
    # corr += sum((correlation_matrix(psi, "S+", "S-")))
    # corr = corr - sum(expect(psi, "Sz"))
    S2_mpo = total_spin_squared_mpo(sites)
    @show S2_val = inner(psi', S2_mpo, psi)

    println("S = ", (-1 + sqrt(1 + 4*S2_val)) / 2)
end 
end


# @time begin
#     psi = load_ψ("wf.h5")

#     corr = 0
#     corr += sum((correlation_matrix(psi, "Sz", "Sz")))
#     corr += sum(0.5*(correlation_matrix(psi, "S+", "S-")))
#     corr += sum(0.5*(correlation_matrix(psi, "S-", "S+")))
#     println(corr)
# end


# println(S2_val)
# println("⟨S²⟩ = ", real(S2_val))
# print("S =", (-1+sqrt(1+4*real(S2_val)))/2)


# sites = siteinds("Electron", 2)
# state = ["Up", "Dn"]
# psi = productMPS(sites, state)

# os = OpSum()
# os += 0.5, "S+", 1, "S-", 2
# os += 0.5, "S-", 1, "S+", 2

# W = MPO(os, sites)
# println("⟨S₁·S₂⟩ spin-flip part: ", inner(psi', W, psi))

# using ITensors

# # Define 2-site system
# sites = siteinds("Electron", 2)

# # Build the singlet state: (|up, dn⟩ - |dn, up⟩)/√2
# up_dn = productMPS(sites, ["Up", "Dn"])
# dn_up = productMPS(sites, ["Dn", "Up"])
# psi = (up_dn - dn_up) / √2

# # Build the spin-flip part of Si · Sj
# os = OpSum()
# os += 0.5, "S+", 1, "S-", 2
# os += 0.5, "S+", 2, "S-", 1
# os += 0.5, "S+", 1, "S-", 2
# os += 0.5, "S+", 2, "S-", 1

# W = MPO(os, sites)
# println("⟨ψ|S+S- + S-S+|ψ⟩ = ", real(inner(psi', W, psi)))
