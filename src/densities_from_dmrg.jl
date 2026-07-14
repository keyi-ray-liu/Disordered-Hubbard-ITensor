#!/usr/bin/env julia
using ITensors, ITensorMPS
using HDF5
using Printf

# -----------------------------
# Settings
# -----------------------------
ROOT    = get(ENV, "ROOT", pwd())
psi_tag = get(ENV, "PSI_TAG", "psi")

phis = (0.24, 0.25, 0.26)

# If true: exclude border (interior only). If false: output all sites.
INTERIOR_ONLY = get(ENV, "INTERIOR_ONLY", "0") == "1"

# Output file
out_file = joinpath(ROOT, "ntot_phi024_025_026_U8_20e.dat")

# -----------------------------
# Helpers
# -----------------------------
flux_folder(φ::Real) = @sprintf("FM_flux%.2f_U8_10_10", φ)
workdir(φ::Real) = joinpath(ROOT, flux_folder(φ), "work")

function load_ψ(work::AbstractString; wf_name="temp_wf1.h5", tag="psi1")
  wf_path = joinpath(work, wf_name)
  isfile(wf_path) || error("Missing wavefunction file:\n  $wf_path")
  ψ = nothing
  h5open(wf_path, "r") do f
    ψ = read(f, tag, MPS)
  end
  return ψ::MPS
end

# Your requested fast path
charge_density(ψ::MPS) = expect(ψ, "Ntot")

# row/col mapping you gave:
#   col = (j - 1) % L(sys) + 1
#   row = div(j - 1, L(sys)) + 1
function colrow(j::Int, Lx::Int)
  col = (j - 1) % Lx + 1
  row = (j - 1) ÷ Lx + 1
  return col, row
end
function infer_LxLy(N::Int)
  if haskey(ENV, "LX") && haskey(ENV, "LY")
    Lx = parse(Int, ENV["LX"])
    Ly = parse(Int, ENV["LY"])
    Lx*Ly == N || error("LX*LY != Nsites: LX=$Lx LY=$Ly -> $(Lx*Ly) but Nsites=$N")
    return Lx, Ly
  end
  s = round(Int, sqrt(N))
  s*s == N || error("Cannot infer LX/LY from N=$N. Set LX and LY env vars.")
  return s, s
end

function is_interior(j::Int, Lx::Int, Ly::Int)
  col, row = colrow(j, Lx)
  return (2 <= col <= Lx-1) && (2 <= row <= Ly-1)
end

# -----------------------------
# Main
# -----------------------------
function main()
  println("ROOT = ", ROOT)
  println("psi_tag = ", psi_tag)
  println("INTERIOR_ONLY = ", INTERIOR_ONLY)

  # Load densities at each phi
  dens = Dict{Float64, Vector{Float64}}()
  Nsites = nothing

  for φ in phis
    wdir = workdir(φ)
    isdir(wdir) || error("Missing work directory:\n  $wdir")

    println("Loading ψ at φ=$(φ) from $wdir ...")
    ψ = load_ψ(wdir; tag=psi_tag)

    println("Computing n = expect(ψ, \"Ntot\") at φ=$(φ) ...")
    n = charge_density(ψ)
    dens[φ] = collect(Float64.(n))
    if Nsites === nothing
      Nsites = length(n)
    elseif length(n) != Nsites
      error("Nsites mismatch across φ: got $(length(n)) vs $Nsites")
    end
  end

  Lx, Ly = infer_LxLy(Nsites)
  println("Using Lx=$Lx, Ly=$Ly (your row/col convention)")

  # Write formatted output
  open(out_file, "w") do io
    @printf(io, "# Site-resolved Ntot from DMRG wf.h5\n")
    @printf(io, "# ROOT=%s\n", ROOT)
    @printf(io, "# psi_tag=%s\n", psi_tag)
    @printf(io, "# Lx=%d Ly=%d\n", Lx, Ly)
    @printf(io, "# Columns: j  col  row  n(phi=%.2f)  n(phi=%.2f)  n(phi=%.2f)\n",
            phis[1], phis[2], phis[3])

    for j in 1:Nsites
      if INTERIOR_ONLY && !is_interior(j, Lx, Ly)
        continue
      end
      col, row = colrow(j, Lx)
      n1 = dens[phis[1]][j]
      n2 = dens[phis[2]][j]
      n3 = dens[phis[3]][j]
      @printf(io, "%4d %3d %3d  % .16e  % .16e  % .16e\n", j, col, row, n1, n2, n3)
    end
  end

  println("Wrote: ", out_file)
end

main()
