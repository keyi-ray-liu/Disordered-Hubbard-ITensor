# !/usr/bin/env julia

using ITensors
using HDF5
using Printf
using Statistics

# -----------------------------
# User settings
# -----------------------------

# Parent directory containing FM_flux_*. If not set, defaults to pwd().
ROOT = get(ENV, "ROOT", pwd())

# HDF5 dataset tag for the MPS inside wf.h5
psi_tag = "psi1"


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


function make_S2_mpo(sites)
    N = length(sites)
    os = OpSum()

    for i in 1:N
        for j in 1:N
            os += 1.0, "Sz", i, "Sz", j
            os += 0.5, "S+", i, "S-", j
            os += 0.5, "S-", i, "S+", j
        end
    end

    return MPO(os, sites)
end

function spin_diagnostics(psi::MPS, sites; cutoff=1e-12, maxdim=2000)
    S2 = make_S2_mpo(sites)

    # <psi|S^2|psi>
    # S2_exp = real(inner(psi, S2, psi))
    S2_exp = real(inner(psi', S2, psi))

    # Apply S^2 to psi
    S2psi = apply(S2, psi; cutoff=cutoff, maxdim=maxdim)
    normalize!(S2psi)

    # Be careful: normalizing S2psi loses the norm information.
    # So instead, for <(S^2)^2>, use inner(S2, psi, S2, psi) if available.
    # S2_sq_exp = real(inner(S2, psi, S2, psi))

    S2psi = apply(S2, psi; cutoff=cutoff, maxdim=maxdim)
    S2_sq_exp = real(inner(S2psi, S2psi))

    variance = S2_sq_exp - S2_exp^2

    # Numerical cleanup
    if abs(variance) < 1e-10
        variance = 0.0
    end

    S_eff = (-1 + sqrt(1 + 4*S2_exp)) / 2

    return (
        S2 = S2_exp,
        S_eff = S_eff,
        S2_squared = S2_sq_exp,
        variance = variance,
    )
end


N = 16
# sites = siteinds("Electron", N; conserve_qns=true)


psi = load_ψ_from_workdir("/Users/ynl42/Documents/GitHub/Disordered-Hubbard-ITensor/src/work/"; wf_name="temp_wf2.h5", tag="psi" )
sites = siteinds(psi)


# Suppose psi is your converged DMRG MPS
diag = spin_diagnostics(psi, sites)

println("<S^2>       = ", diag.S2)
println("S_eff       = ", diag.S_eff)
println("<(S^2)^2>   = ", diag.S2_squared)
println("Var(S^2)    = ", diag.variance)