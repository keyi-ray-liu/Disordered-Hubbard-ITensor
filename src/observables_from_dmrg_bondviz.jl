#!/usr/bin/env julia
# observables_from_dmrg_bondviz.jl
#
# Pure-Julia plotting (Plots.jl/GR), with:
#   - Densities: colored circles at sites
#   - Correlators: colored bonds between nearest-neighbor sites (no direction)
#   - Currents: quiver arrows (optional), as before
#
# Dependencies: ITensors, HDF5, Statistics, Plots
#
# One-time install:
#   using Pkg; Pkg.add(["ITensors","HDF5","Plots"])
#
# Example run (from directory containing wf.h5):
#   julia observables_from_dmrg_bondviz.jl
#
# Example run (from anywhere):
#   julia -e 'include("observables_from_dmrg_bondviz.jl");
#             run_observables(wf_path="/wrk/ynl42/Bfield_7x7/FM_flux_0.25_U_1/work/wf.h5",
#                             Lx=7, Ly=7, flux=0.25, outdir="figs", savefigs=true)'

using ITensors, ITensorMPS
using HDF5
using Statistics
using Printf
using Plots

gr()

# -------------------------
# IO: load MPS from wf.h5
# -------------------------
function load_ψ(wf_path::AbstractString; tag::AbstractString="psi1")
  isfile(wf_path) || error("Missing wf.h5 file: $wf_path")
  h5open(wf_path, "r") do f
    try
      return read(f, tag, MPS)
    catch err
      keys_ = try collect(keys(f)) catch; String[] end
      msg = "Failed to read MPS tag=\"$tag\" from:\n  $wf_path\n"
      if !isempty(keys_) 
       msg *= "Top-level keys in HDF5 file: $(join(keys_, ", "))\n"
      end
      msg *= "Original error: $(err)"
      error(msg)
    end
  end
end

# -------------------------
# Lattice helpers (your convention)
# -------------------------
coord(n, Lx) = ((n-1) % Lx + 1, (n-1) ÷ Lx + 1)  # (x,y) with y=row (1-based)

# nearest-neighbor bonds (open boundaries)
function nn_bonds(Lx::Int, Ly::Int)
  N = Lx * Ly
  bonds = Tuple{Int,Int}[]
  for i in 1:N
    xi, yi = coord(i, Lx)
    if xi < Lx
      push!(bonds, (i, i+1))
    end
    if yi < Ly
      push!(bonds, (i, i+Lx))
    end
  end
  return bonds
end

# -------------------------
# Densities
# -------------------------
charge_density(ψ::MPS) = expect(ψ, "Ntot")
spin_density_sz(ψ::MPS) = expect(ψ, "Sz")

# -------------------------
# Connected correlator helper
# -------------------------
function connected(C::AbstractMatrix, a::AbstractVector, b::AbstractVector=a)
  @assert size(C,1) == length(a)
  @assert size(C,2) == length(b)
  return C .- (a * b')
end

# -------------------------
# Print correlator helper
# -------------------------
function print_nn_bonds(C, Lx, Ly; label="", imag_tol=1e-12)
  bonds = nn_bonds(Lx, Ly)
  println("\nNearest-neighbor $label:")
  for (i,j) in bonds
    v = C[i,j]
    if abs(imag(v)) > imag_tol
      @printf("(%d,%d)–(%d,%d): Re= %+ .6e  Im= %+ .2e  |Im|>tol\n",
              coord(i,Lx)..., coord(j,Lx)..., real(v), imag(v))
    else
      @printf("(%d,%d)–(%d,%d): %+ .6e\n",
              coord(i,Lx)..., coord(j,Lx)..., real(v))
    end
  end
end


# -------------------------
# Correlators (Electron op names)
# -------------------------
corr_SzSz(ψ::MPS) = correlation_matrix(ψ, "Sz", "Sz")

# Spin–spin correlators ⟨Sᵢ·Sⱼ⟩, ⟨SᵢᶻSⱼᶻ⟩
# Physical meaning of the sign: 
# + → ferromagnetic 
# - → antiferromagnetic 

function corr_SdotS(ψ::MPS)
  SzSz = correlation_matrix(ψ, "Sz", "Sz")
  SpSm = correlation_matrix(ψ, "S+", "S-")
  SmSp = correlation_matrix(ψ, "S-", "S+")
  return SzSz .+ 0.5 .* (SpSm .+ SmSp)
end

# Connected spin correlators ⟨Sᵢ·Sⱼ⟩₍c₎
# Physical meaning: “Is there correlation beyond mean-field?”
# + → enhanced correlation
# - → suppressed correlation
# 0 → trivial product state behavior

function corr_SdotS_connected(ψ::MPS, SdotS::AbstractMatrix)
  sz = expect(ψ, "Sz")
  sp = expect(ψ, "S+")
  sm = expect(ψ, "S-")
  prod = (sz * sz') .+ 0.5 .* (sp * sm' .+ sm * sp')
  return SdotS .- prod
end

# Density–density correlators ⟨nᵢnⱼ⟩ and ⟨nᵢnⱼ⟩₍c₎
# ⟨nᵢnⱼ⟩: reflects filling and Hartree physics
# connected ⟨nᵢnⱼ⟩₍c₎ : 
  # + → clustering / attraction
  # - → anticorrelation / repulsion
corr_NN(ψ::MPS) = correlation_matrix(ψ, "Ntot", "Ntot")

# -------------------------
# Currents (your gauge + formula)
# -------------------------
function peierls_phase(i, j, Lx, B)
  xi, yi = coord(i, Lx)
  xj, yj = coord(j, Lx)
  if yi == yj && abs(xj - xi) == 1
    row = yj
    return -2π * B * row
  end
  return 0.0
end

function bond_current(i, j, tij, C_up, C_dn)
  cup_ij = C_up[i,j]; cup_ji = C_up[j,i]
  cdn_ij = C_dn[i,j]; cdn_ji = C_dn[j,i]
  cij = cup_ij + cdn_ij
  cji = cup_ji + cdn_ji
  return im * (tij * cij - conj(tij) * cji)
end

function compute_bond_currents(ψ; Lx, Ly, t=-1.0, flux=0.0)
  N = Lx * Ly
  println("Computing correlation matrices for currents…")
  C_up = correlation_matrix(ψ, "Cdagup", "Cup")
  C_dn = correlation_matrix(ψ, "Cdagdn", "Cdn")

  Jdict = Dict{Tuple{Int,Int}, ComplexF64}()

  for i in 1:N
    xi, yi = coord(i, Lx)

    if xi < Lx
      j = i + 1
      ϕ = peierls_phase(i, j, Lx, flux)
      tij = t * exp(-1im * ϕ)
      Jdict[(i,j)] = bond_current(i, j, tij, C_up, C_dn)
    end

    if yi < Ly
      j = i + Lx
      ϕ = peierls_phase(i, j, Lx, flux)
      tij = t * exp(-1im * ϕ)
      Jdict[(i,j)] = bond_current(i, j, tij, C_up, C_dn)
    end
  end

  return Jdict
end

"""
    site_correlator(C; part=:real, symavg=false)

Compute per-site averaged correlator:
  c[i] = (1/(N-1)) * sum_{j≠i} C[i,j]

Options:
- part=:real (default) uses real(C[i,j]); :abs uses abs(C[i,j]); :full uses C[i,j] (complex output).
- symavg=true averages C[i,j] with C[j,i] before summing (helps suppress tiny anti-Hermitian noise).
"""
function site_correlator(C::AbstractMatrix; part::Symbol=:real, symavg::Bool=false)
  N = size(C, 1)
  @assert size(C, 2) == N

  out = (part == :full) ? zeros(eltype(C), N) : zeros(Float64, N)

  for i in 1:N
    s = zero(eltype(C))
    for j in 1:N
      j == i && continue
      v = symavg ? 0.5*(C[i,j] + C[j,i]) : C[i,j]
      s += v
    end

    s /= (N - 1)

    if part == :full
      out[i] = s
    elseif part == :abs
      out[i] = abs(s)
    else
      out[i] = real(s)
    end
  end

  return out
end

# -------------------------
# Plotting: sites as colored circles
# -------------------------
function save_site_circles(vals::AbstractVector, Lx::Int, Ly::Int;
                           title::AbstractString, outfile::AbstractString,
                           ms::Real=10, cmap=:Purples, diverging=false)

  N = Lx * Ly
  @assert length(vals) == N

  xs = Float64[]; ys = Float64[]; zs = Float64[]
  for i in 1:N
    x,y = coord(i, Lx)
    push!(xs, x); push!(ys, y); push!(zs, real(vals[i]))
  end

  zsf = filter(isfinite, zs) # drops NaN and Inf
  clims = diverging ? (-maximum(abs.(zsf)), maximum(abs.(zsf))) :
                      (minimum(zsf), maximum(zsf))

  p = scatter(xs, ys;
              marker_z=zs,
              ms=ms,
              msw=0,
              xlabel="x", ylabel="y",
              title=title,
              xlim=(0, Lx+1), ylim=(0, Ly+1),
              aspect_ratio=:equal,
              legend=false,
              color=cmap,
              clims=clims,
              colorbar=true)

  savefig(p, outfile)
  return p
end

# -------------------------
# Plotting: correlators as colored bonds (nearest neighbors)
# -------------------------
function save_bond_colormap(C::AbstractMatrix, Lx::Int, Ly::Int;
                            title::AbstractString, outfile::AbstractString,
                            lw::Real=4, symavg::Bool=true, cmap=:coolwarm, 
                            diverging=true)

  N = Lx * Ly
  @assert size(C,1) == N && size(C,2) == N

  bonds = nn_bonds(Lx, Ly)

  xs = Float64[]; ys = Float64[]; zs = Float64[]
  for (i,j) in bonds
    xi, yi = coord(i, Lx)
    xj, yj = coord(j, Lx)

    v = symavg ? 0.5*(real(C[i,j]) + real(C[j,i])) : real(C[i,j])

    push!(xs, xi); push!(ys, yi); push!(zs, v)
    push!(xs, xj); push!(ys, yj); push!(zs, v)
    push!(xs, NaN); push!(ys, NaN); push!(zs, NaN)
  end


  zsf = filter(isfinite, zs) # drops NaN and Inf
  clims = diverging ? (-maximum(abs.(zsf)), maximum(abs.(zsf))) :
                      (minimum(zsf), maximum(zsf))

  p = plot(xs, ys;
           line_z=zs,
           linewidth=lw,
           xlabel="x", ylabel="y",
           title=title,
           xlim=(0, Lx+1), ylim=(0, Ly+1),
           aspect_ratio=:equal,
           legend=false,
           color=cmap,
           clims=clims,
           colorbar=true)

  # overlay faint lattice points for reference
  xs0 = Float64[]; ys0 = Float64[]
  for i in 1:N
    x,y = coord(i, Lx)
    push!(xs0, x); push!(ys0, y)
  end
  scatter!(p, xs0, ys0; ms=3, msw=0, alpha=0.6)

  savefig(p, outfile)
  return p
end

# -------------------------
# Plotting: currents as quiver (arrow color ∝ |J|)
# -------------------------

function save_current_arrows_color(Jdict, Lx, Ly; outfile="currents_color.png", L0=1, lw=3, cmap=:Purples)
    # Build a single polyline with NaN separators, and a matching line_z vector.
    xs = Float64[]; ys = Float64[]; zs = Float64[]

    for ((i,j), J) in Jdict
        xi, yi = coord(i, Lx)
        xj, yj = coord(j, Lx)

        # midpoint
        xm = (xi + xj)/2
        ym = (yi + yj)/2

        # direction based on sign(real(J))
        if real(J) ≥ 0
            dx = xj - xi
            dy = yj - yi
        else
            dx = xi - xj
            dy = yi - yj
        end

        # normalize
        dnorm = sqrt(dx^2 + dy^2)
        dx /= dnorm
        dy /= dnorm

        # fixed arrow length
        L = L0 *0.47

        x2 = xm + dx*L
        y2 = ym + dy*L
        c  = abs(J)

        # segment (xm,ym) -> (x2,y2), then NaN separator
        push!(xs, xm); push!(ys, ym); push!(zs, c)
        push!(xs, x2); push!(ys, y2); push!(zs, c)
        push!(xs, NaN); push!(ys, NaN); push!(zs, NaN)
    end

    # p = plot(; title="Bond currents",
    #           xlim=(0, Lx+1), ylim=(0, Ly+1),
    #           aspect_ratio=:equal, xlabel="x", ylabel="y", legend=false)

    # # Arrowheads for line segments + colormap from line_z
    # plot!(p, xs, ys;
    #       line_z=zs,
    #       c=cmap,
    #       colorbar=true,
    #       linewidth=lw,
    #       arrow=:arrow)

    # savefig(p, outfile)

    p = plot(; 
    title="Bond currents",
    xlim=(0, Lx+1), ylim=(0, Ly+1),
    aspect_ratio=:equal, xlabel="x", ylabel="y",
    legend=false, grid=false, framestyle=:none,
    background_color=:transparent,
    background_color_inside=:white,
    background_color_subplot=:white
    )
    zmax = maximum(zs)
    clims = zmax > 0 ? (0, zmax) : (0, 1e-12)

    # Arrowheads for line segments + colormap from line_z
    plot!(p, xs, ys;
        line_z=zs,
        c=cmap,
        colorbar=true,
        linewidth=lw,
        arrow=:arrow,
        clims = clims)

    savefig(p, outfile)

    return p
end


# -------------------------
# Driver
# -------------------------
function run_observables(;
    wf_path::AbstractString = "wf.h5",
    psi_tag::AbstractString = "psi1",
    Lx::Int = 7,
    Ly::Int = 7,
    t_hop::Float64 = -1.0,
    flux::Float64 = 0.25,
    outdir::AbstractString = ".",
    savefigs::Bool = true
)
    savefigs && (isdir(outdir) || mkpath(outdir))

    N = Lx * Ly
    println("=== observables_from_dmrg (bond viz) ===")
    println("wf_path = $wf_path")
    println("psi_tag = $psi_tag")
    println("Lx×Ly   = $Lx×$Ly (N=$N)")
    println("t       = $t_hop")
    println("flux    = $flux")
    println("outdir  = $outdir")
    println("savefigs= $savefigs")
    println()

    ψ = load_ψ(wf_path; tag=psi_tag)

    # ---- densities: colored circles ----
    n  = charge_density(ψ)
    sz = spin_density_sz(ψ)

    println("Site observables:")
    for i in 1:length(n)
        x, y = coord(i, Lx)
        @printf("site %3d  (x=%d,y=%d):  n = %.6f   Sz = %.6f\n",
                i, x, y, n[i], sz[i])
    end

    if savefigs
      save_site_circles(n,  Lx, Ly; title="Charge density ⟨Ntot⟩",
                        outfile=joinpath(outdir, "density_charge_sites.png"), ms=12,cmap=:Blues)
      println("Saved: ", joinpath(outdir, "density_charge_sites.png"))

      save_site_circles(sz, Lx, Ly; title="Spin density ⟨Sz⟩",
                        outfile=joinpath(outdir, "density_Sz_sites.png"), ms=12, cmap=:coolwarm, diverging=true)
      println("Saved: ", joinpath(outdir, "density_Sz_sites.png"))

      h5open("densities_temp.h5","w") do f
        f["N"] = n
        f["Sz"] = sz
      end

    end

    # ---- correlators: colored bonds ----
    println("Computing correlators with correlation_matrix…")

    SzSz   = corr_SzSz(ψ)
    SzSz_c = connected(SzSz, sz)

    SdotS   = corr_SdotS(ψ)
    SdotS_c = corr_SdotS_connected(ψ, SdotS)

    # c_SdotS_site   = site_correlator(SzSz;   part=:real, symavg=true)
    c_SdotSC_site  = site_correlator(SzSz_c; part=:real, symavg=true)

    println("\nSite-averaged correlator (S dot S connected):")
    for i in 1:length(c_SdotSC_site)
      x,y = coord(i, Lx)
      @printf("site %3d (x=%d,y=%d):  %.6e\n", i, x, y, c_SdotSC_site[i])
    end

    NN   = corr_NN(ψ)
    NN_c = connected(NN, n)
  
    # print_nn_bonds(SzSz, Lx, Ly; label="SzSz")
    # print_nn_bonds(SzSz_c, Lx, Ly; label="SzSz (connected)")
    # print_nn_bonds(SdotS, Lx, Ly; label="SdotS")
    # print_nn_bonds(SdotS_c, Lx, Ly; label="SdotS (connected)")
    # print_nn_bonds(NN, Lx, Ly; label="NN")
    # print_nn_bonds(NN_c, Lx, Ly; label="NN (connected)")

    if savefigs
      save_bond_colormap(SzSz,   Lx, Ly; title="⟨Sz_i Sz_j⟩ on NN bonds",
                         outfile=joinpath(outdir, "bond_SzSz.png"), lw=5, cmap=:coolwarm)
      println("Saved: ", joinpath(outdir, "bond_SzSz.png"))

      save_bond_colormap(SzSz_c, Lx, Ly; title="⟨Sz_i Sz_j⟩_connected on NN bonds",
                         outfile=joinpath(outdir, "bond_SzSz_connected.png"), lw=5, cmap=:coolwarm)
      println("Saved: ", joinpath(outdir, "bond_SzSz_connected.png"))

      save_bond_colormap(SdotS,   Lx, Ly; title="⟨S_i · S_j⟩ on NN bonds",
                         outfile=joinpath(outdir, "bond_SdotS.png"), lw=5, cmap=:coolwarm)
      println("Saved: ", joinpath(outdir, "bond_SdotS.png"))

      save_bond_colormap(SdotS_c, Lx, Ly; title="⟨S_i · S_j⟩_connected on NN bonds",
                         outfile=joinpath(outdir, "bond_SdotS_connected.png"), lw=5, cmap=:coolwarm)
      println("Saved: ", joinpath(outdir, "bond_SdotS_connected.png"))

      save_bond_colormap(NN,   Lx, Ly; title="⟨N_i N_j⟩ on NN bonds",
                         outfile=joinpath(outdir, "bond_NN.png"), lw=5, cmap=:Blues, diverging=false)
      println("Saved: ", joinpath(outdir, "bond_NN.png"))

      save_bond_colormap(NN_c, Lx, Ly; title="⟨N_i N_j⟩_connected on NN bonds",
                         outfile=joinpath(outdir, "bond_NN_connected.png"), lw=5, cmap=:Blues, diverging=false)
      println("Saved: ", joinpath(outdir, "bond_NN_connected.png"))

      h5open("correlators_temp.h5","w") do f
        f["Czz"] = SzSz
        f["Cdot"] = SdotS
        f["Czz_c"] = SzSz_c
        f["Cdot_c"] = SdotS_c
        f["NN"] = NN
        f["NN_c"] = NN_c
      end
      # # to load later 
      # using HDF5
      # h5open("correlators.h5","r") do f
      # Czz = read(f["Czz"])
      # end
    end

    # ---- currents (arrows) ----
    println("Computing bond currents…")
    Jdict = compute_bond_currents(ψ; Lx=Lx, Ly=Ly, t=t_hop, flux=flux)
    println("Computed $(length(Jdict)) bond currents.")
    println("Bond currents (sorted):")
    for ((i,j), J) in sort(collect(Jdict); by=x->x[1])
        @printf("J(%3d -> %3d) = % .6e  +  i% .6e   |J|=% .6e\n",
                i, j, real(J), imag(J), abs(J))
    end


    if savefigs
      save_current_arrows_color(Jdict, Lx, Ly; outfile=joinpath(outdir, "currents_quiver.png"), L0=1,cmap=:Purples)
      println("Saved: ", joinpath(outdir, "currents_quiver.png"))
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
  run_observables()
end
