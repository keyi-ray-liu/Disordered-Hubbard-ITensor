using DelimitedFiles
using Glob
using ITensors
using ITensors.HDF5
using ITensors: OneITensor, check_hascommoninds, linkind, siteinds, tr
using JSON3
using LinearAlgebra
using StatsBase
using Suppressor
using ITensorTDVP
using Random

include("systems.jl")
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
include("LSR_SIAM.jl")
include("QEsystems.jl")
include("NF.jl")

disable_blas = true

if disable_blas
  BLAS.set_num_threads(1)
  ITensors.Strided.set_num_threads(1)
end 

if ARGS != []
    test = false
else
    test = true
end 

if test

    #DPT_wrapper()
    #NF_wrapper()
    QE_SIAM_wrapper()


else

    ARG = ARGS[1]

    if ARG == "DPT"
        DPT_wrapper()

    elseif ARG == "occ"
        dyna_occ()

    elseif ARG == "EE" 
        dyna_EE()

    elseif ARG == "dptcurrent"
        dyna_dptcurrent()

    elseif ARG == "NF_suqre"
        NF_wrapper()

    elseif ARG == "QE_SIAM"
        QE_SIAM_wrapper()

    else
        error("Unrecognized option")
    end 
end 
#run_Static_chain(12, 6)
#run_QE_two(1.0, 12, 6, true; staticex = 5, dp=1.0, init="Both")
#run_DPT(2.0, 2, 2, 1.0, 2.0)
#DPT_wrapper()
#LSR_SIAM_wrapper()

