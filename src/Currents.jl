using DelimitedFiles
using Glob
using ITensors
using ITensors: OneITensor, linkind, siteinds, tr
using JSON3
using LinearAlgebra
using StatsBase
using Suppressor
using ITensorTDVP
using Random
using HDF5
using Logging

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


filename = "wf.h5"
using HDF5
# print(HDF5.ishdf5(filename))
ψ = load_ψ(filename)

# dens = expect(ψ, "Ntot")
# # print(dens)
# updens, dndens = expect(ψ, "Nup", "Ndn")
# # print(updens, dndens)

# zzcorr = correlation_matrix(ψ,"Sz","Sz") # Sx or Sy dosen't work since it does not conserve e- number

# zzcorr = correlation_matrix(ψ,"I","I")
# print(zzcorr[1, :])
# for i in 1:9
#     print(real(zzcorr[i, :]), ' ')
# end




"""
  weighted_corr(ψ, sites; t, φ)

Compute the matrix T[i,j] = ⟨ψ|t[i,j]*exp(1im*φ(i,j))*c†up(i)*cup(j)|ψ⟩
for an MPS ψ on fermion sites `sites`.  
- t can be either a scalar or an N×N array of magnitudes.  
- φ(i,j) should return a real-valued phase (Float64).
"""

function CurrentMagnetic(j::Int, N::Int, L::Int, t, flux::Float64)

    hop = []

    # column and row number 
    col = (j - 1) % L + 1

    # flux = 0.1
    # calculate Peierls phase, default is Landau Gauge, only exist on vertical bond
    phase = exp(1im * 2 * π * flux * col )

    # not at end of col
    if j % L != 0
        append!(hop, [[t, t, j + 1]])
    end 
    
    # not at end of row
    if j <= N - L
        append!(hop, [[t*phase, t*phase, j + L]])
    end 
    
    # println(hop)
    return hop
end


function landaux(xi::Number, yi::Number, xj::Number, yj::Number, flux::Number)
    y_avg   = (yi + yj) / 2
    delta_x = xj - xi
    phase   = -2 * π * flux * y_avg * delta_x
    return exp(1im * phase)
end

function landauy(xi::Number, yi::Number, xj::Number, yj::Number, flux::Number)
    x_avg   = (xi + xj) / 2
    delta_y = yj - yi
    phase   = 2 * π * flux * x_avg * delta_y
    return exp(1im * phase)
end

function weighted_corr_test(ψ::MPS; t=-1.0, L=3, N=9, flux=0.1)
    corr_up = correlation_matrix(ψ,"Cup","Cdagup")
    corr_dn = correlation_matrix(ψ,"Cdn","Cdagdn")

    corr_total = corr_up + corr_dn
    current_x, current_y = [], []
    for j in 1:N
        col = (j-1)%L + 1
        row = div(j-1, L) + 1
    
        println("which j = ", j)
        # current x
        if col == 1 # only to right 
            # phase = landauy(col, row, col+1, row, flux)
            phase = landaux(col, row, col+1, row, flux)
            Jx = 1im*(t*phase*corr_total[j, j+1] - conj(t*phase)*corr_total[j+1, j])
        elseif col == L # only to left
            # phase = landauy(col, row, col-1, row, flux)
            phase = landaux(col, row, col-1, row, flux)
            Jx = 1im*(t*phase*corr_total[j, j-1] - conj(t*phase)*corr_total[j-1, j])
        else # both left and right 
            # phaser = landauy(col, row, col+1, row, flux)
            phaser = landaux(col, row, col+1, row, flux)
            Jxr = 1im*(t*phaser*corr_total[j, j+1] - conj(t*phaser)*corr_total[j+1, j])

            # phasel = landauy(col, row, col-1, row, flux)
            phasel = landaux(col, row, col-1, row, flux)
            Jxl = 1im*(t*phasel*corr_total[j, j-1] - conj(t*phasel)*corr_total[j-1, j])
            Jx = Jxr - Jxl
        end
        append!(current_x, Jx)

        # current y
        if row == 1 # only to up 
            # phase = landauy(col, row, col, row+1, flux)
            phase = landaux(col, row, col, row+1, flux)
            Jy = 1im*(t*phase*corr_total[j, j+L] - conj(t*phase)*corr_total[j+L, j])
        elseif row == L # only to down
            # phase = landauy(col, row, col, row-1, flux)
            phase = landaux(col, row, col, row-1, flux)
            Jy = 1im*(t*phase*corr_total[j, j-L] - conj(t*phase)*corr_total[j-L, j])
        else # both up and down 
            # phaser = landauy(col, row, col, row+1, flux)
            phaser = landaux(col, row, col, row+1, flux)
            Jyu = 1im*(t*phaser*corr_total[j, j+L] - conj(t*phaser)*corr_total[j+L, j])

            # phasel = landauy(col, row, col, row-1, flux)
            phasel = landaux(col, row, col, row-1, flux)
            Jyd = 1im*(t*phasel*corr_total[j, j-L] - conj(t*phasel)*corr_total[j-L, j])
            Jy = Jyu - Jyd
        end
        append!(current_y, Jy)
    end
    return real(current_x), real(current_y)
end
    
function weighted_corr(ψ::MPS, sites; t=1.0, L=3, N=9, flux=0.1)

    T = zeros(ComplexF64, N, N)
    operators = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]

    for j in 1:N
        for (v..., k) in CurrentMagnetic(j, N, L, t, flux)
            k = trunc(Int, real(k)) # hop to pair 
            # println("j, k, v = ", j, k, v)
            for (i, operator) in enumerate(operators)
                # println(v[i])
                if v[i] != 0
                    
                    op1, op2 = operator

                    res = OpSum()
                    # println("one")
                    # opsum += (v[i],"Cdagup",4,"Cup",5)
                    # println("aaa")
                    res += v[i], op1, j, op2,  k
                    # println("1111")

                    res += conj(v[i]), op1, k, op2, j
                    # println("operator", res)

                    O = MPO(res, siteinds(ψ))
                    T[j,k] = inner(ψ', O, ψ)
                    # ex_W = inner(psi',W,psi)
                    # println("two")

                end 
            end 

        end 
    end 
    return T

end


# println("test first")
# total_num = 9
sites = siteinds("Electron", 9) # N = Lx*Ly
# println("test second")

tij = -1.0
# L = 3
# flux = 0.1        

# println("test third")

# 4) compute the weighted correlator
# T = weighted_corr(ψ, sites; t=tij, L=3, N=9, flux=0.1)
# println(T[1,:])
Jx, Jy = weighted_corr_test(ψ; t=-1.0, L=3, N=9, flux=0.1)
println(Jx)
println(Jy)
# total = up + dn
# println(total[1,:])
JX =[-0.10012813996613934, -0.20025627993227746, 0.1001281399661381, 
3.3133218391157016e-16, 4.597017211338539e-16, -1.2836953722228372e-16,
 0.1001281399661401, 0.20025627993227818, -0.10012813996613808]
JY =[0.10012813996614034, 2.0469737016526324e-16, -0.10012813996613856, 
0.20025627993228018, 1.3964523981613297e-16, -0.20025627993227602, 
-0.10012813996613984, 6.505213034913027e-17, 0.10012813996613744]

# current_up = correlation_matrix(ψ, "Cdagup", "Cup")
# current_dn = correlation_matrix(ψ, "Cdagdn", "Cdn")
# for i in 1:9

#     println("up = ", current_up[i,:])
#     println("dn = ", current_dn[i,:])
#     println(" ", )
# end