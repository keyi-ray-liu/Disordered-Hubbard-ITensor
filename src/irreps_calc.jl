# !/usr/bin/env julia

using ITensors
using HDF5
using Printf
using Statistics
using ITensorMPS


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



# load the file 
# N = 16
# HDF5 dataset tag for the MPS inside wf.h5
psi_tag = "psi"
wf_tag = "temp_wf1.h5"
psi = load_ψ_from_workdir("/Users/ynl42/Documents/GitHub/Disordered-Hubbard-ITensor/src/work/"; wf_name=wf_tag, tag=psi_tag)
sites = siteinds(psi)

## calculate the total spin directly 

let
    corr = 0
    corr += sum((correlation_matrix(psi, "Sz", "Sz")))
    corr += sum(0.5*(correlation_matrix(psi, "S+", "S-")))
    corr += sum(0.5*(correlation_matrix(psi, "S-", "S+")))
    println("Total S2 = ", corr)
    println("Total S = ", (-1 + sqrt(1 + 4 * corr)) / 2)
end

## functions to calculate irreps
function ITensors.op(::OpName"FSWAP", ::SiteType"Electron",
                     s1::Index, s2::Index)

    G = ITensor(prime(s1), prime(s2), dag(s1), dag(s2))

    # ITensor Electron basis:
    # 1 = Emp
    # 2 = Up
    # 3 = Dn
    # 4 = UpDn
    parity = [0, 1, 1, 0]

    for a in 1:dim(s1)
        for b in 1:dim(s2)
            sign = isodd(parity[a] * parity[b]) ? -1.0 : 1.0

            # Input:  a on site 1, b on site 2
            # Output: b on site 1, a on site 2
            G[prime(s1) => b,
              prime(s2) => a,
              dag(s1) => a,
              dag(s2) => b] = sign
        end
    end

    return G
end

# coordinates helper 
function site_index(x, y, L)
    return x + L*y + 1
end

function coord(i, L)
    x = (i - 1) % L
    y = (i - 1) ÷ L
    return x, y
end

function make_perm(L, mapxy)
    N = L * L
    p = zeros(Int, N)

    for i in 1:N
        x, y = coord(i, L)
        xp, yp = mapxy(x, y)
        p[i] = site_index(xp, yp, L)
    end

    return p
end

function d4_perms(L)
    return Dict(
        "E"     => collect(1:L*L),

        # C4 counterclockwise: (x,y) -> (L-1-y,x)
        "C4"    => make_perm(L, (x,y) -> (L - 1 - y, x)),

        # C2: (x,y) -> (L-1-x,L-1-y)
        "C2"    => make_perm(L, (x,y) -> (L - 1 - x, L - 1 - y)),

        # C4^{-1}: (x,y) -> (y,L-1-x)
        "C4m"   => make_perm(L, (x,y) -> (y, L - 1 - x)),

        # vertical reflection: x -> L-1-x
        "sig_v" => make_perm(L, (x,y) -> (L - 1 - x, y)),

        # horizontal reflection: y -> L-1-y
        "sig_h" => make_perm(L, (x,y) -> (x, L - 1 - y)),

        # main diagonal reflection: (x,y) -> (y,x)
        "sig_d" => make_perm(L, (x,y) -> (y, x)),

        # anti-diagonal reflection: (x,y) -> (L-1-y,L-1-x)
        "sig_a" => make_perm(L, (x,y) -> (L - 1 - y, L - 1 - x))
    )
end

function adjacent_swaps_for_perm(p)
    # p[old_site] = new_site.
    # After applying U_g, final position new_site contains old_site g^{-1}(new_site).
    target = invperm(p)

    current = collect(1:length(p))
    swaps = Int[]

    for pos in 1:length(p)
        j = findfirst(==(target[pos]), current)

        while j > pos
            push!(swaps, j - 1)
            current[j - 1], current[j] = current[j], current[j - 1]
            j -= 1
        end
    end

    return swaps
end

function apply_site_perm_fswap(psi::MPS, sites, p;
                               cutoff=1e-10,
                               maxdim=2000)

    ψg = copy(psi)

    swaps = adjacent_swaps_for_perm(p)

    for b in swaps
        G = op("FSWAP", sites[b], sites[b + 1])
        ψg = apply(G, ψg; cutoff=cutoff, maxdim=maxdim)
        ψg = noprime(ψg)
    end

    return ψg
end

function symmetry_overlap_fswap(psi::MPS, sites, p;
                                cutoff=1e-10,
                                maxdim=2000)

    ψg = apply_site_perm_fswap(psi, sites, p;
                               cutoff=cutoff,
                               maxdim=maxdim)

    return inner(psi, ψg)
end


# #### example for C4
# L = 4
# perms = d4_perms(L)

# o_C4 = symmetry_overlap_fswap(psi, sites, perms["C4"];
#                               cutoff=1e-9,
#                               maxdim=2000)

# println("<psi|C4|psi> = ", o_C4)

# for all D4 
function d4_overlaps(psi::MPS, sites, L;
                     cutoff=1e-10,
                     maxdim=2000)

    perms = d4_perms(L)

    overlaps = Dict{String,ComplexF64}()

    for name in ["E", "C4", "C2", "C4m", "sig_v", "sig_h", "sig_d", "sig_a"]
        if name == "E"
            overlaps[name] = inner(psi, psi)
        else
            overlaps[name] = symmetry_overlap_fswap(psi, sites, perms[name];
                                                    cutoff=cutoff,
                                                    maxdim=maxdim)
        end
    end

    return overlaps
end

# overlaps = d4_overlaps(psi, sites, L;
#                        cutoff=1e-8,
#                        maxdim=1000)

# for name in ["E", "C4", "C2", "C4m", "sig_v", "sig_h", "sig_d", "sig_a"]
#     println(name, "  ", overlaps[name])
# end

function d4_irrep_weights_from_overlaps(o)
    e    = real(o["E"])
    c4   = real(o["C4"])
    c2   = real(o["C2"])
    c4m  = real(o["C4m"])
    sv   = real(o["sig_v"])
    sh   = real(o["sig_h"])
    sd   = real(o["sig_d"])
    sa   = real(o["sig_a"])

    weights = Dict{String,Float64}()

    weights["A1"] = (e + c4 + c2 + c4m + sv + sh + sd + sa) / 8

    weights["A2"] = (e + c4 + c2 + c4m - sv - sh - sd - sa) / 8

    weights["B1"] = (e - c4 + c2 - c4m + sv + sh - sd - sa) / 8

    weights["B2"] = (e - c4 + c2 - c4m - sv - sh + sd + sa) / 8

    weights["E"]  = 0.5 * (e - c2)

    return weights
end

L = 4
overlaps = d4_overlaps(psi, sites, L;
                       cutoff=1e-8,
                       maxdim=1000)

weights = d4_irrep_weights_from_overlaps(overlaps)

println("D4 overlaps:")
for k in ["E", "C4", "C2", "C4m", "sig_v", "sig_h", "sig_d", "sig_a"]
    println(k, " = ", overlaps[k])
end

println("\nD4 irrep weights:")
for ir in ["A1", "A2", "B1", "B2", "E"]
    println(ir, " = ", weights[ir])
end