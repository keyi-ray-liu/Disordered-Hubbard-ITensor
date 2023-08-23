# load variables specifically in linsolve
using ITensors: linkind, siteinds, check_hascommoninds, OneITensor
#using ITensors: @debug_check, @timeit_debug, @printf
#using KrylovKit: eigsolve, linsolve
#import ITensors: permute, position!
using ITensorTDVP
# now load the rest of the functions 
using ITensors
using DelimitedFiles
using Suppressor, Glob
using ITensors.HDF5
using LinearAlgebra
using ITensorGaussianMPS: correlation_matrix_to_mps, slater_determinant_to_mps

include("main.jl")
include("setpara.jl")
include("setdisorder.jl")
include("utils.jl")
include("init_state.jl")
include("init_ham.jl")
include("search.jl")
include("correlation.jl")
include("data_gen.jl")
include("process.jl")
include("shift_and_invert.jl")
include("timeevolve.jl")
include("ham_helper.jl")
include("wrappers.jl")
#include("projection.jl")

"""The top level function that controls workflow. Use operation mode to select code function"""
function top()

  # test handle for more streamlined testing environment in julia REPL
  #test = false
  test = false
  disable_blas = true

  if disable_blas
    BLAS.set_num_threads(1)
    ITensors.Strided.set_num_threads(1)
  end 

  if test

    println("TEST TEST TEST")

    paras = setpara()
    main(paras;)
    #eigensolver_wrapper()
    #QEdyna_wrapper()
    #QE("2", "0.0855")
    #eigen_overlap()
    #NF("0.01", "1", "40", "10", "5", "chain")
    #QE_dynamic()
    #paras = setpara(L=12, N=6, CN=6, ex=3, int_ee=2.0, int_ne=2.0, guess=false, method="DMRG", sweepdim=100, 
    #sweepcnt=40, noise=false, QE=2, QN=true, QEen=0.6, dp= [-1.0, 1.0] , ζ_dp = [0.5, 0.5] , QEloc = [[-2.0], [13.0]])
    #sweepcnt=40, noise=false, QE=0, QN=true, QEen=0.0, dp= [] , ζ_dp = [] , QEloc = [])
    #main(paras)


    
  else
    if ARGS == []
      println("must specify operating mode.")

    elseif ARGS[1] == "0"
      plasmonstat( ARGS[2], ARGS[3])

    elseif ARGS[1] == "1"
      cal_observe()

    elseif ARGS[1] == "2"
      GS_wrapper()

    elseif ARGS[1] == "3"
      truedisorder( ARGS[2], ARGS[3] )

    # test NF
    elseif ARGS[1] == "4"
      NF(ARGS[2], ARGS[3], ARGS[4], ARGS[5], ARGS[6], ARGS[7])

    # testing QE
    elseif ARGS[1] == "5"
      QE_wrapper(ARGS[2], ARGS[3])
    
    elseif ARGS[1] == "6"
      QEdyna_wrapper()

    elseif ARGS[1] == "7"
      temp_occ(ARGS[2])

    elseif ARGS[1] == "8"
      cal_overlap()

    elseif ARGS[1] == "9"
      eigen_overlap()

    elseif ARGS[1] == "10"
      eigensolver_wrapper()

    elseif ARGS[1] == "11"
      eigenplot()

    elseif ARGS[1] == "12"
      time_corr_plot()
      
    else
      println("not a valid operating mode")
    
    end 
  end 
end 

top()

