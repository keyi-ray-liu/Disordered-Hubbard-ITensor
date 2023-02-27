# load variables specifically in linsolve
using ITensors: linkind, siteinds, check_hascommoninds, OneITensor
#using ITensors: @debug_check, @timeit_debug, @printf
#using KrylovKit: eigsolve, linsolve
#import ITensors: permute, position!
using ITensorTDVP
# now load the rest of the functions 
using ITensors
using DelimitedFiles
using Suppressor
using ITensors.HDF5
using LinearAlgebra
using ITensorGaussianMPS

include("main.jl")
include("setpara.jl")
include("setdisorder.jl")
include("utils.jl")
include("init_state.jl")
include("init_ham.jl")
include("search.jl")
include("correlation.jl")
include("operations.jl")
include("shift_and_invert.jl")
include("timeevolve.jl")
#include("projection.jl")

"""The top level function that controls workflow. Use operation mode to select code function"""
function top()

  # test handle for more streamlined testing environment in julia REPL
  #test = false
  test = true
  disable_blas = true

  if disable_blas
    BLAS.set_num_threads(1)
    ITensors.Strided.set_num_threads(1)
  end 

  if test

    println("TEST TEST TEST")
    QE_dynamic()
    #paras = setpara(L=60, N=30, CN=30, ex=20, int_ee=2.0, int_ne=2.0, guess=false, method="DMRG", sweepdim=300, 
    #sweepcnt=50, noise=false, QE=2, QN=true, QEen=0.6, dp= 1.0 * [1.0, -1.0], ζ_dp = [0.5, 0.5], QEoffset=1.0)
    #main(paras)
    
  else
    if ARGS == []
      println("must specify operating mode. 0: plasmon (L, lam), 1: gscc, 2: scan stat, 3: true disorder stat, 4: single instance of main, 5: QE, 6: time-evolution")

    elseif ARGS[1] == "0"
      plasmonstat( ARGS[2], ARGS[3])

    elseif ARGS[1] == "1"
      cal_observe()

    elseif ARGS[1] == "2"
      scandisorder( ARGS[2], ARGS[3])

    elseif ARGS[1] == "3"
      truedisorder( ARGS[2], ARGS[3] )

    elseif ARGS[1] == "4"
      paras = setpara(L=[3, 30], N=45, ex=20, int_ee=1.0, guess=false, sweepdim=2000, sweepcnt=100, noise=false)
      @time main(paras)

    # testing QE
    elseif ARGS[1] == "5"
      QE(ARGS[2], ARGS[3])
    
    elseif ARGS[1] == "6"
      QE_dynamic()

    else
      println("not a valid operating mode")
    
    end 
  end 
end 

top()

