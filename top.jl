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
#include("projection.jl")

"""The top level function that controls workflow. Use operation mode to select code function"""
function top()

  # test handle for more streamlined testing environment in julia REPL
  test = true

  if test
    main(L=12, N=6, ex=20, int_ee=1.0, guess=false, method="DMRG", sweepdim=200, sweepcnt=200, noise=false, QE=0, QN=false)

  else
    if ARGS == []
      println("must specify operating mode. 0: plasmon (L, lam), 1: gscc, 2: scan stat, 3: true disorder stat, 4: single instance of main")

    elseif ARGS[1] == "0"
      plasmonstat( ARGS[2], ARGS[3])

    elseif ARGS[1] == "1"
      cal_observe()

    elseif ARGS[1] == "2"
      scandisorder( ARGS[2], ARGS[3])

    elseif ARGS[1] == "3"
      truedisorder( ARGS[2], ARGS[3] )

    elseif ARGS[1] == "4"
      @time main(L=[3, 30], N=45, ex=20, int_ee=1.0, guess=false, sweepdim=2000, sweepcnt=100, noise=false)

    else
      println("not a valid operating mode")
    
    end 
  end 
end 

top()

