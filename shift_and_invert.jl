
#function IndexSet_ignore_missing(is::Union{Index,Nothing}...)
#  return IndexSet(filter(i -> i isa Index, is))
#end

"""
    shift_and_invert(H::MPO,ϕ::MPS, ψ0::MPS;kwargs...)

                    
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
function shift_and_invert(H::MPO, ϕ::MPS,  λ::Float64, sites, para::Dict)
  
  L = para["L"]
  N = para["N"]

  # generate a random guess for the new state
  # This is the opposite of the init_state function initialization
  state = append!([ "Emp" for n=1:N] , ["Occ" for n=1:L -N])
  ψ0 = randomMPS(sites,state)

  ψ0 = ϕ
  return linsolve(H, ϕ, ψ0, -λ, 1.0)

end


