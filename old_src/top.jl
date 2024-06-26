# load variables specifically in linsolve
using DelimitedFiles
using Glob
using ITensorGaussianMPS: correlation_matrix_to_mps, slater_determinant_to_mps
using ITensorTDVP
using ITensors
using ITensors.HDF5
using ITensors: OneITensor, check_hascommoninds, linkind, siteinds, tr
using JSON3
using LinearAlgebra
using StatsBase
using Suppressor



const QESITES = 2

#include("behavior_test.jl")
include("main.jl")
include("setpara.jl")
include("setdisorder.jl")
include("utils.jl")
include("init_state.jl")
include("init_ham.jl")
include("search.jl")
include("observable.jl")
include("data_gen.jl")
include("process.jl")
include("shift_and_invert.jl")
include("timeevolve.jl")
include("ham_hop.jl")
include("ham_coul.jl")
include("ham_onsite.jl")
include("ham_qe.jl")
include("ham_qn.jl")
include("ham_sd.jl")
include("wrappers.jl")
include("QE.jl")
include("SD_close.jl")
include("io.jl")
include("time_obs.jl")
include("gqs.jl")
#include("projection.jl")

"""The top level function that controls workflow. Use operation mode to select code function"""
function top()

  # test handle for more streamlined testing environment in julia REPL
  #test = false
  
  if ARGS != []
    test = false
  else
    test = true
  end 

  disable_blas = true

  if disable_blas
    BLAS.set_num_threads(1)
    ITensors.Strided.set_num_threads(1)
  end 

  if test

    println("TEST TEST TEST")
    
    #rm("work/", recursive=true)
    #GQS_dyna_wrapper()
    #benchmark()
    #iter_sd_wrapper()
    #source_drain_wrapper()
    #transport_wrapper()
    #paras = setpara(;L=6, N=3, ex=6, QE=0, dynamode="left", headoverride=0); main(paras)
    #eigensolver_wrapper()
    #time_obs_wrapper( "occ")
    #QEdyna_wrapper()
    #eigen_overlap()
    NF_wrapper()
    #time_obs_wrapper("GQS")
    #REPL_test_wrapper()
    #corr_wrapper()
    #QE_wrapper("2", "1.0")

    
  else

    if ARGS[1] == "0"
      plasmonstat( ARGS[2], ARGS[3])

    elseif ARGS[1] == "1"
      cal_observe()

    elseif ARGS[1] == "2"
      GS_wrapper()

    elseif ARGS[1] == "3"
      truedisorder( ARGS[2], ARGS[3] )

    # test NF
    elseif ARGS[1] == "4"
      NF_wrapper()

    # testing QE
    elseif ARGS[1] == "5"
      QE_wrapper(ARGS[2], ARGS[3])
    
    elseif ARGS[1] == "6"
      QEdyna_wrapper()

    elseif ARGS[1] == "7"
      time_obs_wrapper( ARGS[2])

    elseif ARGS[1] == "8"
      cal_overlap()

    elseif ARGS[1] == "9"
      eigen_overlap()

    elseif ARGS[1] == "10"
      eigensolver_wrapper()

    elseif ARGS[1] == "11"
      eigenplot()

    elseif ARGS[1] == "12"
      corr_wrapper()

    elseif ARGS[1] == "13"
      source_drain_wrapper()

    elseif ARGS[1] == "14"
      time_tcd()

    elseif ARGS[1] == "15"
      transport_wrapper()

    elseif ARGS[1] == "16"
      GQS_dyna_wrapper()

    else
      println("not a valid operating mode")
    
    end 
  end 
end 

top()

