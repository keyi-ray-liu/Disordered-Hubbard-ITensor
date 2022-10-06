
#function IndexSet_ignore_missing(is::Union{Index,Nothing}...)
#  return IndexSet(filter(i -> i isa Index, is))
#end

function permute(
  M::AbstractMPS, ::Tuple{typeof(linkind),typeof(siteinds),typeof(linkind)}
)::typeof(M)
  M̃ = typeof(M)(length(M))
  for n in 1:length(M)
    lₙ₋₁ = linkind(M, n - 1)
    lₙ = linkind(M, n)
    s⃗ₙ = sort(Tuple(siteinds(M, n)); by=plev)
    M̃[n] = permute(M[n], filter(!isnothing, (lₙ₋₁, s⃗ₙ..., lₙ)))
  end
  set_ortho_lims!(M̃, ortho_lims(M))
  return M̃
end

"""
    shift_and_invert(H::MPO,H2::MPO,psi0::MPS;kwargs...)
    shift_and_invert(H::MPO,H2::MPO,psi0::MPS,sweeps::Sweeps;kwargs...)
                    
Use the density matrix renormalization group (shift_and_invert) algorithm
to optimize a matrix product state (MPS) such that it is the
eigenvector of lowest eigenvalue of a Hermitian matrix H,
represented as a matrix product operator (MPO).
The MPS `psi0` is used to initialize the MPS to be optimized.

The number of sweeps of thd shift_and_invert algorithm is controlled by
passing the `nsweeps` keyword argument. The keyword arguments
`maxdim`, `cutoff`, `noise`, and `mindim` can also be passed
to control the cost versus accuracy of the algorithm - see below
for details.

Alternatively the number of sweeps and accuracy parameters can
be passed through a `Sweeps` object, though this interface is 
no longer preferred.

Returns:
* `psi::MPS` - optimized MPS

Keyword arguments:
* `nsweeps::Int` - number of "sweeps" of shift_and_invert to perform

Optional keyword arguments:
* `maxdim` - integer or array of integers specifying the maximum size allowed for the bond dimension or rank of the MPS being optimized
* `cutoff` - float or array of floats specifying the truncation error cutoff or threshold to use for truncating the bond dimension or rank of the MPS
* `noise` - float or array of floats specifying strength of the "noise term" to use to aid convergence
* `mindim` - integer or array of integers specifying the minimum size of the bond dimension or rank, if possible
* `outputlevel::Int = 1` - larger outputlevel values make shift_and_invert print more information and 0 means no output
* `observer` - object implementing the [Observer](@ref observer) interface which can perform measurements and stop shift_and_invert early
* `write_when_maxdim_exceeds::Int` - when the allowed maxdim exceeds this value, begin saving tensors to disk to free memory in large calculations
"""
function shift_and_invert(H::MPO, H2::MPO, psi0::MPS, sweeps::Sweeps, lambda::Float64; kwargs...)
  check_hascommoninds(siteinds, H, psi0)
  check_hascommoninds(siteinds, H, psi0')
  check_hascommoninds(siteinds, H2, psi0)
  check_hascommoninds(siteinds, H2, psi0')
  # Permute the indices to have a better memory layout
  # and minimize permutations
  H = permute(H, (linkind, siteinds, linkind))
  PH = ProjMPO(H)
  PH2 = ProjMPO(H2)
  return shift_and_invert(PH, PH2, psi0, sweeps, lambda; kwargs...)
end





function shift_and_invert(PH, PH2, psi0::MPS, sweeps::Sweeps, lambda::Float64; kwargs...)
  if length(psi0) == 1
    error(
      "`shift_and_invert` currently does not support system sizes of 1. You can diagonalize the MPO tensor directly with tools like `LinearAlgebra.eigen`, `KrylovKit.linsolve`, etc.",
    )
  end

  @debug_check begin
    # Debug level checks
    # Enable with ITensors.enable_debug_checks()
    checkflux(psi0)
    checkflux(PH)
  end

  which_decomp::Union{String,Nothing} = get(kwargs, :which_decomp, nothing)
  svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")
  obs = get(kwargs, :observer, NoObserver())
  outputlevel::Int = get(kwargs, :outputlevel, 1)

  write_when_maxdim_exceeds::Union{Int,Nothing} = get(
    kwargs, :write_when_maxdim_exceeds, nothing
  )

  # linsolve kwargs
  linsolve_tol::Number = get(kwargs, :linsolve_tol, 1e-14)
  linsolve_krylovdim::Int = get(kwargs, :linsolve_krylovdim, 3)
  linsolve_maxiter::Int = get(kwargs, :linsolve_maxiter, 5)
  linsolve_verbosity::Int = get(kwargs, :linsolve_verbosity, 0)

  # TODO: add support for non-Hermitian shift_and_invert
  ishermitian::Bool = get(kwargs, :ishermitian, true)

  # TODO: add support for targeting other states with shift_and_invert
  # (such as the state with the largest eigenvalue)
  # get(kwargs, :linsolve_which_eigenvalue, :SR)

  # TODO: use this as preferred syntax for passing arguments
  # to linsolve
  #default_linsolve_args = (tol = 1e-14, krylovdim = 3, maxiter = 1,
  #                         verbosity = 0, ishermitian = true,
  #                         which_eigenvalue = :SR)
  #linsolve = get(kwargs, :linsolve, default_linsolve_args)

  # Keyword argument deprecations
  if haskey(kwargs, :maxiter)
    error("""maxiter keyword has been replaced by linsolve_krylovdim.
             Note: compared to the C++ version of ITensor,
             setting linsolve_krylovdim 3 is the same as setting
             a maxiter of 2.""")
  end

  if haskey(kwargs, :errgoal)
    error("errgoal keyword has been replaced by linsolve_tol.")
  end

  if haskey(kwargs, :quiet)
    error("quiet keyword has been replaced by outputlevel")
  end

  psi = copy(psi0)

  # reference state phi0, unchanged through out
  phi0 = copy(psi0)

  N = length(psi)

  if !isortho(psi) || orthocenter(psi) != 1
    orthogonalize!(psi, 1)
  end
  @assert isortho(psi) && orthocenter(psi) == 1

  # we set the unprojected site to 1 on the HdagH term, i.e. LHS
  position!(PH2, psi, 1)

  # we set the unprojected site to 1 on the H term, i.e. RHS
  position!(PH, psi, phi0, 1)

  for sw in 1:nsweep(sweeps)
    
    sw_time = @elapsed begin
      maxtruncerr = 0.0

      if !isnothing(write_when_maxdim_exceeds) &&
        maxdim(sweeps, sw) > write_when_maxdim_exceeds
        if outputlevel >= 2
          println(
            "write_when_maxdim_exceeds = $write_when_maxdim_exceeds and maxdim(sweeps, sw) = $(maxdim(sweeps, sw)), writing environment tensors to disk",
          )
        end
        PH2 = disk(PH2)
      end

      for (b, ha) in sweepnext(N)
        @debug_check begin
          checkflux(psi)
          checkflux(PH)
          checkflux(PH2)
        end

        # set both RHS and LHS to site b
        @timeit_debug timer "shift_and_invert: position!" begin
          RHS = position!(PH, psi, phi0, b)
          position!(PH2, psi, b)
        end

        @debug_check begin
          checkflux(psi)
          checkflux(PH)
          checkflux(PH2)
        end
        

        @timeit_debug timer "shift_and_invert: psi[b]*psi[b+1]" begin
          guess = psi[b] * psi[b + 1]
        end
        
        # for linsolve, LHS = PH2, RHS = PH. Initial guess still provided as phi (not to be confused with phi0, which is unchanged)
        @timeit_debug timer "shift_and_invert: linsolve" begin
          vecs, conv = linsolve(
            PH2,
            RHS,
            -lambda, 1;
            ishermitian=ishermitian,
            tol=linsolve_tol,
            krylovdim=linsolve_krylovdim,
            maxiter=linsolve_maxiter,
          )
          
        end
        
        print("linsolve passed")

        #print(typeof(vecs))
        #print("vec", vecs, "\n\n\n", "conv", conv)
        
        # These next steps in replacing the old state should be the same 
        newpsi::ITensor = vecs

        ortho = ha == 1 ? "left" : "right"

        drho = nothing
        if noise(sweeps, sw) > 0.0
          @timeit_debug timer "shift_and_invert: noiseterm" begin
            # Use noise term when determining new MPS basis
            drho = noise(sweeps, sw) * noiseterm(PH, newpsi, ortho)
          end
        end

        @debug_check begin
          checkflux(newpsi)
        end

        @timeit_debug timer "shift_and_invert: replacebond!" begin
          spec = replacebond!(
            psi,
            b,
            newpsi;
            maxdim=maxdim(sweeps, sw),
            mindim=mindim(sweeps, sw),
            cutoff=cutoff(sweeps, sw),
            eigen_perturbation=drho,
            ortho=ortho,
            normalize=true,
            which_decomp=which_decomp,
            svd_alg=svd_alg,
          )
        end
        maxtruncerr = max(maxtruncerr, spec.truncerr)

        @debug_check begin
          checkflux(psi)
          checkflux(PH)
        end

        if outputlevel >= 2
          @printf("Sweep %d, half %d, bond (%d,%d)\n", sw, ha, b, b + 1)
          @printf(
            "  Truncated using cutoff=%.1E maxdim=%d mindim=%d\n",
            cutoff(sweeps, sw),
            maxdim(sweeps, sw),
            mindim(sweeps, sw)
          )
          @printf(
            "  Trunc. err=%.2E, bond dimension %d\n", spec.truncerr, dim(linkind(psi, b))
          )
          flush(stdout)
        end

        sweep_is_done = (b == 1 && ha == 2)

        linsolvemeasure!(
          obs;
          psi=psi,
          bond=b,
          sweep=sw,
          half_sweep=ha,
          spec=spec,
          outputlevel=outputlevel,
          sweep_is_done=sweep_is_done,
        )
      end
    end
    if outputlevel >= 1
      @printf(
        "After sweep %d  maxlinkdim=%d maxerr=%.2E time=%.3f\n",
        sw,
        maxlinkdim(psi),
        maxtruncerr,
        sw_time
      )
      flush(stdout)
    end
    isdone = linsolvecheckdone!(obs; psi=psi, sweep=sw, outputlevel=outputlevel)
    isdone && break
  end
  return psi
end

function _shift_and_invert_sweeps(;
  nsweeps, maxdim=typemax(Int), mindim=1, cutoff=1E-8, noise=0.0, kwargs...
)
  sweeps = Sweeps(nsweeps)
  setmaxdim!(sweeps, maxdim...)
  setmindim!(sweeps, mindim...)
  setcutoff!(sweeps, cutoff...)
  setnoise!(sweeps, noise...)
  return sweeps
end

function shift_and_invert(x1, x2, psi0::MPS; kwargs...)
  return shift_and_invert(x1, x2, psi0, _shift_and_invert_sweeps(; kwargs...); kwargs...)
end

function shift_and_invert(x1, psi0::MPS; kwargs...)
  return shift_and_invert(x1, psi0, _shift_and_invert_sweeps(; kwargs...); kwargs...)
end
