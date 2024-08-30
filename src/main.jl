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

#using DataGraphs: edge_data, vertex_data
#using Dictionaries: Dictionary
# using Graphs: nv, vertices, edges, src, dst
# #using ITensors: ITensors
# #using ITensors.ITensorMPS: ITensorMPS
# using ITensorNetworks:
#   ITensorNetworks,
#   OpSum,
#   ttn,
#   apply,
#   dmrg,
#   inner,
#   linkdims,
#   mpo,
#   random_mps,
#   random_ttn,
#   siteinds

#using ITensorNetworks : tdvp
#using ITensorNetworks.ModelHamiltonians: ModelHamiltonians
#using KrylovKit: eigsolve
# using NamedGraphs:  rem_vertex!, add_vertex!, add_edge!, NamedGraph
# using Observers: observer
using Test: @test, @test_broken, @testset
using ITensorUnicodePlots: @visualize
using StableRNGs: StableRNG


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


const DISABLE_BLAS = true

if DISABLE_BLAS
  BLAS.set_num_threads(1)
  ITensors.Strided.set_num_threads(1)
end 

if ARGS != []
    test = false
else
    test = true
end 

if test
    map(rm, Glob.glob( "corr*", getworkdir()))
    #rm(getworkdir(), force=true, recursive=true)
    #plot_mix()
    #static_tcd()
    #chain_wrapper()
    #prepare_wavepacket()
    #test_embedding_wrapper()
    QE_gaussian_wrapper()
    #gen_GS_scan()
    #dyna_pϕ()
    #GQS_wrapper()
    #set_SD()
    #DPT_wrapper()
    #SD_wrapper()
    #solve_QE()
    #DPT_corr()
    #benchmark()
    #dyna_dptcurrent_mix()
    #NF_wrapper()
    #dyna_EE()
    #QE_SIAM_wrapper()
    #QE_wrapper("QE_two")
    #QE_HOM_wrapper()
    #DPT_graph_test()
    #runtest()
    #solve_QE()
    #dyna_tcd()

    return nothing


else

    ARG = ARGS[1]

    if ARG == "DPT"
        DPT_wrapper()

    elseif ARG == "occ"
        dyna_occ()

    elseif ARG == "EE" 
        dyna_EE()

    elseif ARG == "DPT_corr"
        DPT_corr()

    elseif ARG == "dyna_TCD"
        dyna_tcd()

    elseif ARG == "mixcurrent"
        dyna_dptcurrent_mix()

    elseif ARG == "dptcurrent"
        dyna_dptcurrent()

    elseif ARG == "NF_square"
        NF_wrapper()

    # elseif ARG == "QE_SIAM"
    #     QE_SIAM_wrapper()


    elseif ARG == "QE_two"
        QE_wrapper("QE_two")

    elseif ARG == "QE_parallel"
        QE_wrapper("QE_HOM")

    elseif ARG == "QE_gaussian"
        QE_gaussian_wrapper()

    elseif ARG == "Solve_QE"
        solve_QE()

    elseif ARG == "plot_mix"
        plot_mix()

    elseif ARG == "TCD"
        static_tcd()

    elseif ARG == "SD"
        SD_wrapper()

    elseif ARG == "chain"
        chain_wrapper()

    else
        error("Unrecognized option")
    end 
end 
#run_Static_chain(12, 6)
#run_QE_two(1.0, 12, 6, true; staticex = 5, dp=1.0, init="Both")
#run_DPT(2.0, 2, 2, 1.0, 2.0)
#DPT_wrapper()
#LSR_SIAM_wrapper()

