using ITensors
using HDF5
using LinearAlgebra
using Printf

# -----------------------------
# Load MPS
# -----------------------------
function load_ψ(wf_path::AbstractString; tag::AbstractString="psi")
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

# -----------------------------
# Geometry helpers
# row-major: i = x + Lx*y + 1
# x,y are 0-based
# -----------------------------
site_index(x, y, Lx) = x + Lx*y + 1

function coords(i, Lx)
    i0 = i - 1
    x = i0 % Lx
    y = div(i0, Lx)
    return x, y
end

# perm[i] = image of site i under symmetry operation g
function c4v_permutations(Lx, Ly)
    Lx == Ly || error("C4v requires Lx == Ly. For rectangular lattices use C2v instead.")
    L = Lx
    N = Lx * Ly

    perms = Dict{String, Vector{Int}}()

    # Identity
    perms["E"] = collect(1:N)

    # C4: (x,y) -> (L-1-y, x)
    perms["C4"] = [site_index(L-1-y, x, Lx) for i in 1:N for (x,y) in (coords(i,Lx),)]

    # C2: (x,y) -> (L-1-x, L-1-y)
    perms["C2"] = [site_index(L-1-x, L-1-y, Lx) for i in 1:N for (x,y) in (coords(i,Lx),)]

    # C4 inverse: (x,y) -> (y, L-1-x)
    perms["C4i"] = [site_index(y, L-1-x, Lx) for i in 1:N for (x,y) in (coords(i,Lx),)]

    # vertical mirror: x -> L-1-x
    perms["sv"] = [site_index(L-1-x, y, Lx) for i in 1:N for (x,y) in (coords(i,Lx),)]

    # horizontal mirror: y -> L-1-y
    perms["sh"] = [site_index(x, L-1-y, Lx) for i in 1:N for (x,y) in (coords(i,Lx),)]

    # diagonal mirror: (x,y) -> (y,x)
    perms["sd"] = [site_index(y, x, Lx) for i in 1:N for (x,y) in (coords(i,Lx),)]

    # anti-diagonal mirror: (x,y) -> (L-1-y, L-1-x)
    perms["sa"] = [site_index(L-1-y, L-1-x, Lx) for i in 1:N for (x,y) in (coords(i,Lx),)]

    return perms
end

# -----------------------------
# Overlap <phi | U(perm) | psi>
#
# Convention:
# perm[i] = image of old site i.
# The ket physical index at old site i is identified with
# the bra physical index at new site perm[i].
# -----------------------------
function overlap_permuted(phi::MPS, psi::MPS, perm::Vector{Int})
    N = length(psi)
    N == length(phi) || error("MPS lengths do not match.")
    length(perm) == N || error("Permutation length does not match MPS length.")

    s_phi = siteinds(phi)
    s_psi = siteinds(psi)

    T = ITensor(1.0)

    # Bra network
    for i in 1:N
        T *= dag(phi[i])
    end

    # Ket network with crossed physical-site contractions
    for i in 1:N
        A = psi[i]
        A = replaceind(A, s_psi[i], s_phi[perm[i]])
        T *= A
    end

    return scalar(T)
end

# -----------------------------
# Build D(g) matrix in degenerate subspace
# D[m,n] = <psi_m | U(g) | psi_n>
# -----------------------------
function symmetry_matrix(states::Vector{MPS}, perm::Vector{Int})
    n = length(states)
    D = zeros(ComplexF64, n, n)

    for m in 1:n
        for n2 in 1:n
            D[m,n2] = overlap_permuted(states[m], states[n2], perm)
        end
    end

    return D
end

# ordinary Gram matrix <psi_m|psi_n>
function gram_matrix(states::Vector{MPS})
    n = length(states)
    S = zeros(ComplexF64, n, n)

    for m in 1:n
        for n2 in 1:n
            S[m,n2] = inner(states[m], states[n2])
        end
    end

    return S
end

# -----------------------------
# C4v character table for individual group elements
# -----------------------------
function c4v_characters()
    return Dict(
        "A1" => Dict("E"=>1, "C4"=>1,  "C2"=>1,  "C4i"=>1,  "sv"=>1,  "sh"=>1,  "sd"=>1,  "sa"=>1),
        "A2" => Dict("E"=>1, "C4"=>1,  "C2"=>1,  "C4i"=>1,  "sv"=>-1, "sh"=>-1, "sd"=>-1, "sa"=>-1),
        "B1" => Dict("E"=>1, "C4"=>-1, "C2"=>1,  "C4i"=>-1, "sv"=>1,  "sh"=>1,  "sd"=>-1, "sa"=>-1),
        "B2" => Dict("E"=>1, "C4"=>-1, "C2"=>1,  "C4i"=>-1, "sv"=>-1, "sh"=>-1, "sd"=>1,  "sa"=>1),
        "E"  => Dict("E"=>2, "C4"=>0,  "C2"=>-2, "C4i"=>0,  "sv"=>0,  "sh"=>0,  "sd"=>0,  "sa"=>0)
    )
end

# -----------------------------
# Main classification
# -----------------------------
function classify_c4v_irrep(wf_paths::Vector{String}; Lx::Int, Ly::Int, tag::String="psi")
    states = [load_ψ(p; tag=tag) for p in wf_paths]

    println("Loaded $(length(states)) states.")

    perms = c4v_permutations(Lx, Ly)
    chars_table = c4v_characters()

    S = gram_matrix(states)
    println("\nGram matrix <psi_m|psi_n>:")
    display(S)

    println("\nDeviation from orthonormality ||S-I|| = ", norm(S - I))

    # If states are not exactly orthonormal, use generalized character:
    # chi(g) = tr(S^{-1} D(g))
    Sinv = inv(S)

    chars_subspace = Dict{String, ComplexF64}()

    println("\nCharacters of your ground-state subspace:")
    for opname in ["E", "C4", "C2", "C4i", "sv", "sh", "sd", "sa"]
        D = symmetry_matrix(states, perms[opname])
        χ = tr(Sinv * D)
        chars_subspace[opname] = χ
        @printf("%4s : % .8f %+.8fi\n", opname, real(χ), imag(χ))
    end

    println("\nIrrep multiplicities from character projection:")
    group_order = 8

    multiplicities = Dict{String, ComplexF64}()

    for irrep in ["A1", "A2", "B1", "B2", "E"]
        mult = 0.0 + 0.0im
        for opname in keys(perms)
            mult += conj(chars_table[irrep][opname]) * chars_subspace[opname]
        end
        mult /= group_order
        multiplicities[irrep] = mult
        @printf("%2s : % .8f %+.8fi\n", irrep, real(mult), imag(mult))
    end

    println("\nInterpretation:")
    println("Multiplicity near 1 means that irrep appears once in your ground-state manifold.")
    println("Multiplicity near 0 means it is absent.")
    println("For example, E near 1 means a twofold E doublet.")
    println("A1 + B1 both near 1 means two separate one-dimensional irreps are degenerate.")

    return chars_subspace, multiplicities
end

function main()
    work_dir = joinpath(@__DIR__, "work_temp")

    wf_paths = abspath.([
        joinpath(work_dir, "temp_wf1.h5"),
        joinpath(work_dir, "temp_wf2.h5"),
        # joinpath(work_dir, "temp_wf3.h5"),
        # joinpath(work_dir, "temp_wf4.h5"),
    ])

    foreach(p -> println(p, "  exists = ", isfile(p)), wf_paths)

    chars, mults = classify_c4v_irrep(wf_paths; Lx=4, Ly=4, tag="psi")
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
    # wf_paths = abspath.(joinpath(dirname(@__DIR__), "work_temp", f) for f in [
    # "temp_wf1.h5","temp_wf2.h5"])
    # # "temp_wf1.h5","temp_wf2.h5","temp_wf3.h5","temp_wf4.h5"])

    # chars, mults = classify_c4v_irrep(wf_paths; Lx=6, Ly=6, tag="psi")
end