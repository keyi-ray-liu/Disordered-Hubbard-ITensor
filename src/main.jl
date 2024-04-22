using DelimitedFiles
using Glob
using ITensors
using ITensors.HDF5
using ITensors: OneITensor, linkind, siteinds, tr
using JSON3
using LinearAlgebra
using StatsBase
using Suppressor
using ITensorTDVP
using Random

#using DataGraphs: edge_data, vertex_data
#using Dictionaries: Dictionary
using Graphs: nv, vertices, edges, src, dst
#using ITensors: ITensors
#using ITensors.ITensorMPS: ITensorMPS
using ITensorNetworks:
  ITensorNetworks,
  OpSum,
  ttn,
  apply,
  dmrg,
  inner,
  linkdims,
  mpo,
  random_mps,
  random_ttn,
  relabel_sites,
  siteinds
#using ITensorNetworks.ModelHamiltonians: ModelHamiltonians
#using KrylovKit: eigsolve
using NamedGraphs: named_comb_tree, rem_vertex!, add_vertex!, add_edge!, NamedGraph
using Observers: observer
using Test: @test, @test_broken, @testset
using ITensorUnicodePlots: @visualize

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
    #dyna_EE()
    #QE_SIAM_wrapper()
    #QE_two_wrapper()
    QE_parallel_wrapper()


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

    elseif ARG == "NF_square"
        NF_wrapper()

    elseif ARG == "QE_SIAM"
        QE_SIAM_wrapper()


    elseif ARG == "QE_two"
        QE_two_wrapper()

    else
        error("Unrecognized option")
    end 
end 
#run_Static_chain(12, 6)
#run_QE_two(1.0, 12, 6, true; staticex = 5, dp=1.0, init="Both")
#run_DPT(2.0, 2, 2, 1.0, 2.0)
#DPT_wrapper()
#LSR_SIAM_wrapper()

