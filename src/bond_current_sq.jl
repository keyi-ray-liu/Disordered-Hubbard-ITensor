
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

########################################################
# --- Helper: coordinate mapping (1D index → 2D coord)
########################################################
coord(n, Lx) = ((n-1) % Lx + 1, (n-1) ÷ Lx + 1)

###############################################################
# Peierls phase for your gauge:
#   φ = -2π B * row(j) for horizontal bonds
#   φ = 0              for vertical bonds
###############################################################
function peierls_phase(i, j, Lx, B)
    xi, yi = coord(i, Lx)
    xj, yj = coord(j, Lx)

    # horizontal bond?
    if yi == yj && abs(xj - xi) == 1
        row = yj
        return -2π * B * row
    end

    # vertical bond (no phase)
    return 0.0
end


########################################################
# --- Current computation: J = i(t ⟨c†_i c_j⟩ − t* ⟨c†_j c_i⟩ )
########################################################
function bond_current(i, j, tij, C_up, C_dn)
    cup_ij = C_up[i,j]; cup_ji = C_up[j,i]
    cdn_ij = C_dn[i,j]; cdn_ji = C_dn[j,i]

    cij = cup_ij + cdn_ij
    cji = cup_ji + cdn_ji

    # println("the current:", i, "-> ", j, ": ", im * (tij * cij - conj(tij) * cji))
    return im * (tij * cij - conj(tij) * cji)
end

###############################################################
# MAIN FUNCTION:
# Compute all nearest-neighbor bond currents
###############################################################
function compute_bond_currents(ψ; Lx, Ly, t=-1.0, flux=0.0)

    N = Lx * Ly

    println("Computing correlation matrices…")
    C_up = correlation_matrix(ψ, "Cdagup", "Cup")
    C_dn = correlation_matrix(ψ, "Cdagdn", "Cdn")

    Jdict = Dict{Tuple{Int,Int}, ComplexF64}()

    for i in 1:N
        xi, yi = coord(i, Lx)

        # ----------------------------
        # Horizontal neighbor (i → i+1)
        # ----------------------------
        if xi < Lx
            j = i + 1

            ϕ = peierls_phase(i, j, Lx, flux)
            tij = t * exp(-1im * ϕ)
            # println(i, "→", j, ":", tij)
            Jdict[(i,j)] = bond_current(i, j, tij, C_up, C_dn)
        end

        # ----------------------------
        # Vertical neighbor (i → i+Lx)
        # ----------------------------
        if yi < Ly
            j = i + Lx

            ϕ = peierls_phase(i, j, Lx, flux)
            tij = t * exp(-1im * ϕ)
            # print(i, "→", j, ϕ)
            # println(i, "→", j, ":", tij)

            Jdict[(i,j)] = bond_current(i, j, tij, C_up, C_dn)
        end
    end

    return Jdict
end

using PyPlot

coord(n, Lx) = ((n-1) % Lx + 1, (n-1) ÷ Lx + 1)

"""
    plot_current_vectorfield_length(Jdict, Lx, Ly; scale=1.0)

Vector field where the *arrow length* encodes |J|.
"""
function plot_current_vectorfield_length(Jdict, Lx, Ly; scale=1.0)

    fig, ax = subplots(figsize=(6,6))
    ax.set_aspect("equal")

    for ((i,j), J) in Jdict
        xi, yi = coord(i, Lx)
        xj, yj = coord(j, Lx)

        # midpoint of bond
        xm = (xi + xj)/2
        ym = (yi + yj)/2

        # direction of bond
        dx = xj - xi
        dy = yj - yi

        # normalize direction
        norm = sqrt(dx^2 + dy^2)
        dx /= norm
        dy /= norm

        # arrow length = |J|
        L = scale * abs(J)

        ax.arrow(xm, ym, dx*L, dy*L,
                 head_width=0.15, head_length=0.2,
                 length_includes_head=true,
                 color="black", alpha=0.8)
    end

    ax.set_xlim(0, Lx+1)
    ax.set_ylim(0, Ly+1)
    ax.set_title("Bond Currents: Arrow length = |J|")
    fig.tight_layout()
    return fig, ax
end



function plot_current_vectorfield_color(Jdict, Lx, Ly; L0=0.5)

    fig, ax = subplots(figsize=(6,6))
    ax.set_aspect("equal")

    # gather |J| for colormap scaling
    Jvals = [abs(J) for (_,J) in Jdict]
    Jmin, Jmax = minimum(Jvals), maximum(Jvals)
    cmap = plt.cm.Purples
    norm = plt.Normalize(Jmin, Jmax)

    for ((i,j), J) in Jdict

        # ---- Coordinates of sites ----
        xi, yi = coord(i, Lx)
        xj, yj = coord(j, Lx)

        # ---- Choose direction based on sign(J) ----
        if real(J) ≥ 0
            # arrow points i → j
            dx = xj - xi
            dy = yj - yi
            xm = (xi + xj)/2
            ym = (yi + yj)/2
        else
            # arrow points j → i
            dx = xi - xj
            dy = yi - yj
            xm = (xi + xj)/2
            ym = (yi + yj)/2
        end

        # ---- Normalize direction ----
        dnorm = sqrt(dx^2 + dy^2)
        dx /= dnorm
        dy /= dnorm

        # ---- Constant arrow length ----
        L = L0*0.5

        ax.arrow(xm, ym, dx*L, dy*L,
                 head_width=0.15, head_length=0.2,
                 length_includes_head=true,
                 color=cmap(norm(abs(J))), alpha=0.9)
    end

    ax.set_xlim(0, Lx+1)
    ax.set_ylim(0, Ly+1)
    # ax.set_title("Bond Currents")
    plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                 label="|J|")

    fig.tight_layout()
    return fig, ax
end


########################################################
# --- Example usage
########################################################

function run_currents()

    ψ = load_ψ("wf.h5")
    Lx, Ly = 2, 8
    t = -1
    Bflux = 0.5   # your B field

    Jdict = compute_bond_currents(ψ; Lx=Lx, Ly=Ly, t=t, flux=Bflux)

    for ((i,j), J) in Jdict
        println("J($i → $j) = $J")
    end

    println("Computed $(length(Jdict)) bond currents.")

    plot_current_vectorfield_color(Jdict, Lx, Ly)
    PyPlot.show()
end

run_currents()

