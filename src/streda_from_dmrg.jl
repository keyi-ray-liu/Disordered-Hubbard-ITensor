# !/usr/bin/env julia
# streda_from_dmrg.jl
#
# Compute a Středa-type response d n / dφ using already-run DMRG wavefunctions.
#
# Assumes folders exist:
#   FM_flux_0.24_U_10/work/wf.h5
#   FM_flux_0.25_U_10/work/wf.h5
#   FM_flux_0.26_U_10/work/wf.h5
#
# Energies file (optional if CHECK_GAP=false):
#   FM_flux_*/work/wfallenergy   (either HDF5 dataset or plain-text)
#
# Usage:
#   julia streda_from_dmrg.jl
# Or set ROOT externally:
#   ROOT=/path/to/parent julia streda_from_dmrg.jl

using ITensors
using HDF5
using Printf
using Statistics

# -----------------------------
# User settings
# -----------------------------

# Parent directory containing FM_flux_*. If not set, defaults to pwd().
ROOT = get(ENV, "ROOT", pwd())

# Choose which sites you want to treat as "bulk" for n_bulk(φ).
# IMPORTANT: set this to match your lattice indexing convention.
# Examples:
#   bulk_sites = [12]              # single center site
#   bulk_sites = [12,13,17,18]     # 2x2 center on a 4x4 with row-major indexing
bulk_sites = Int[25]  # <-- FILL ME (empty => will use average over ALL sites)

# HDF5 dataset tag for the MPS inside wf.h5
psi_tag = "psi1"

# Flux values to use for the derivative around phi0=0.25
phi0   = 0.25
phi_lo = 0.24
phi_hi = 0.26

# Convert ∂n/∂φ to Chern estimate. Common convention: C ≈ 2π * ∂n/∂φ when φ is flux/flux_quantum per plaquette.
# If your φ is defined differently, adjust this factor.
chern_factor = 1.0

# -----------------------------
# Gap-check control
# -----------------------------
CHECK_GAP = false   # set to false to completely disable gap checking

# -----------------------------
# Robust path helpers
# -----------------------------

flux_folder(φ::Real) = @sprintf("FM_flux_%.2f_U_10", φ)
workdir(φ::Real) = joinpath(ROOT, flux_folder(φ), "work")

function assert_exists(path::AbstractString; kind="path")
  if !ispath(path)
    error("Missing $(kind):\n  $path")
  end
  return path
end

# -----------------------------
# Wavefunction loader (absolute/portable)
# -----------------------------
function load_ψ_from_workdir(work::AbstractString; wf_name::AbstractString="wf.h5", tag::AbstractString="psi1")
  wf_path = joinpath(work, wf_name)
  assert_exists(wf_path; kind="wavefunction file")

  ψ = nothing
  h5open(wf_path, "r") do f
    # ITensors stores MPS in HDF5; `read(f, tag, MPS)` is the standard API.
    try
      ψ = read(f, tag, MPS)
    catch err
      # Helpful diagnostics: list top-level keys if tag is wrong
      keys_ = try
        collect(keys(f))
      catch
        String[]
      end
      msg = "Failed to read MPS tag=\"$tag\" from:\n  $wf_path\n"
      if !isempty(keys_)
      msg *= "Top-level keys in HDF5 file: $(join(keys_, ", "))\n"
      end
      msg *= "Original error: $(err)"
      error(msg)
    end
  end

  return ψ::MPS
end

# -----------------------------
# Energy loader + gap check
# -----------------------------
function read_energies(work::AbstractString; energy_name::AbstractString="wfallenergy")
  ep = joinpath(work, energy_name)
  assert_exists(ep; kind="energy file")

  # Try HDF5 first (common if you saved energies into an HDF5 dataset/file)
  if endswith(lowercase(ep), ".h5") || endswith(lowercase(ep), ".hdf5")
    Es = Float64[]
    h5open(ep, "r") do f
      # Heuristic: if dataset name "E" or "energies" exists, use it; else try first dataset.
      cand = String[]
      for nm in ("E", "energies", "wfallenergy", "energy")
        if haskey(f, nm)
          push!(cand, nm)
        end
      end
      dset = !isempty(cand) ? cand[1] : (length(keys(f)) > 0 ? first(collect(keys(f))) : nothing)
      dset === nothing && error("No datasets found in energy HDF5 file:\n  $ep")
      data = read(f[dset])
      Es = vec(Float64.(data))
    end
    return Es
  end

  # Otherwise treat it as plain-text:
  # Accept formats like:
  #   E0
  #   E0 E1 E2 ...
  #   i  Ei
  #   i  Ei  <other columns>
  nums = Float64[]
  open(ep, "r") do io
    for ln in eachline(io)
      s = strip(ln)
      isempty(s) && continue
      startswith(s, "#") && continue
      parts = split(s)
      for p in parts
        v = tryparse(Float64, p)
        v === nothing && continue
        push!(nums, v)
      end
    end
  end

  if isempty(nums)
    error("Parsed zero numeric values from energy file:\n  $ep\nIf this is HDF5, rename it to *.h5 or adjust read_energies().")
  end

  # Heuristic: if the file is (index, energy, ...) rows, then energies are likely every 2nd token
  if length(nums) ≥ 4
    if abs(nums[1] - 1) < 1e-9 && abs(nums[3] - 2) < 1e-9
      Es = Float64[]
      for k in 2:2:length(nums)
        push!(Es, nums[k])
      end
      return Es
    end
  end

  return sort(nums)
end

function gap_check(work::AbstractString)
  Es = read_energies(work)
  Es_sorted = sort(Es)
  if length(Es_sorted) < 2
    return (gap=NaN, E0=Es_sorted[1], E1=NaN, ok=false)
  end
  E0, E1 = Es_sorted[1], Es_sorted[2]
  gap = E1 - E0
  ok = isfinite(gap) && gap > 0
  return (gap=gap, E0=E0, E1=E1, ok=ok)
end

# -----------------------------
# Density extraction using ITensor utilities
# -----------------------------
function n_profile(ψ::MPS; opname::AbstractString="Ntot")
  # ITensors.jl provides expect(ψ, "OpName") for site-wise expectation values.
  # For Electron sites, common names: "Nup", "Ndn", "Ntot".
  try
    return expect(ψ, opname)
  catch err
    # Try constructing Ntot = Nup + Ndn if needed
    if opname == "Ntot"
      nup = expect(ψ, "Nup")
      ndn = expect(ψ, "Ndn")
      return nup .+ ndn
    else
      rethrow(err)
    end
  end
end

function n_bulk_from_profile(n::AbstractVector{<:Real}, bulk_sites::Vector{Int})
  if isempty(bulk_sites)
    return mean(n)
  else
    @assert all(1 .<= bulk_sites .<= length(n)) "bulk_sites has indices outside 1:$((length(n)))"
    return mean(n[bulk_sites])
  end
end

# -----------------------------
# Core: compute n_bulk(φ), with optional gap check
# -----------------------------
function compute_n_bulk(φ::Real;
                        psi_tag::AbstractString="psi1",
                        bulk_sites::Vector{Int}=Int[],
                        check_gap::Bool=CHECK_GAP)

  wdir = workdir(φ)
  assert_exists(wdir; kind="work directory")

  ψ = load_ψ_from_workdir(wdir; wf_name="wf.h5", tag=psi_tag)
  n = n_profile(ψ; opname="Ntot")
  nb = n_bulk_from_profile(n, bulk_sites)

  if check_gap
    g = gap_check(wdir)
    return (
      n_bulk = nb,
      n_profile = n,
      gap = g.gap,
      E0 = g.E0,
      E1 = g.E1,
      gap_ok = g.ok,
      wdir = wdir
    )
  else
    return (
      n_bulk = nb,
      n_profile = n,
      gap = NaN,
      E0 = NaN,
      E1 = NaN,
      gap_ok = true, # forced true so downstream code doesn't warn
      wdir = wdir
    )
  end
end

# -----------------------------
# Main computation
# -----------------------------
function main()
  @printf("ROOT: %s\n", ROOT)
  @printf("CHECK_GAP: %s\n\n", string(CHECK_GAP))

  res_lo = compute_n_bulk(phi_lo; psi_tag=psi_tag, bulk_sites=bulk_sites)
  res_0  = compute_n_bulk(phi0;   psi_tag=psi_tag, bulk_sites=bulk_sites)
  res_hi = compute_n_bulk(phi_hi; psi_tag=psi_tag, bulk_sites=bulk_sites)

  @printf("Bulk definition: %s\n",
          isempty(bulk_sites) ? "ALL sites (mean density)" : "sites = $(bulk_sites)")

  if CHECK_GAP
    @printf("\nφ = %.2f  n_bulk = %.12f  gap = %.12e  (E0=%.12f, E1=%.12f)  gap_ok=%s\n",
            phi_lo, res_lo.n_bulk, res_lo.gap, res_lo.E0, res_lo.E1, string(res_lo.gap_ok))
    @printf("φ = %.2f  n_bulk = %.12f  gap = %.12e  (E0=%.12f, E1=%.12f)  gap_ok=%s\n",
            phi0,  res_0.n_bulk,  res_0.gap,  res_0.E0,  res_0.E1,  string(res_0.gap_ok))
    @printf("φ = %.2f  n_bulk = %.12f  gap = %.12e  (E0=%.12f, E1=%.12f)  gap_ok=%s\n",
            phi_hi, res_hi.n_bulk, res_hi.gap, res_hi.E0, res_hi.E1, string(res_hi.gap_ok))
  else
    @printf("\nφ = %.2f  n_bulk = %.12f\n", phi_lo, res_lo.n_bulk)
    @printf("φ = %.2f  n_bulk = %.12f\n", phi0,  res_0.n_bulk)
    @printf("φ = %.2f  n_bulk = %.12f\n", phi_hi, res_hi.n_bulk)
  end

  # Central difference around 0.25 using (0.26 - 0.24)/0.02
  dphi = phi_hi - phi_lo
  dn_dphi_central = (res_hi.n_bulk - res_lo.n_bulk) / dphi
  C_est_central = chern_factor * dn_dphi_central

  # One-sided derivatives around 0.25 for a sanity check
  dn_dphi_left  = (res_0.n_bulk  - res_lo.n_bulk) / (phi0 - phi_lo)
  dn_dphi_right = (res_hi.n_bulk - res_0.n_bulk)  / (phi_hi - phi0)
  C_left  = chern_factor * dn_dphi_left
  C_right = chern_factor * dn_dphi_right

  @printf("\n--- Středa response around φ0=%.2f ---\n", phi0)
  @printf("central  dn/dφ ≈ (n(%.2f)-n(%.2f)) / %.3f = %.12e\n", phi_hi, phi_lo, dphi, dn_dphi_central)
  @printf("left     dn/dφ ≈ (n(%.2f)-n(%.2f)) / %.3f = %.12e\n", phi0,  phi_lo, phi0-phi_lo, dn_dphi_left)
  @printf("right    dn/dφ ≈ (n(%.2f)-n(%.2f)) / %.3f = %.12e\n", phi_hi, phi0,  phi_hi-phi0, dn_dphi_right)

  @printf("\nWith chern_factor = %.12f (default 2π):\n", chern_factor)
  @printf("C_est (central) = %.12f\n", C_est_central)
  @printf("C_est (left)    = %.12f\n", C_left)
  @printf("C_est (right)   = %.12f\n", C_right)

  if CHECK_GAP && !(res_lo.gap_ok && res_0.gap_ok && res_hi.gap_ok)
    @printf("\nWARNING: At least one gap check did not pass (gap <= 0 or not finite).\n")
    @printf("If wfallenergy stores only E0 or a nonstandard format, adjust read_energies().\n")
  end

  return nothing
end

main()
