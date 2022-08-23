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

"""The top level function that controls workflow. Use operation mode to select code function"""
function top()

  if ARGS == []
    println("must specify operating mode. 0: plasmon, 1: gscc, 2: scan stat, 3: true disorder stat, 4: single instance of main")

  elseif ARGS[1] == "0"
    plasmonstat( ARGS[2])

  elseif ARGS[1] == "1"
    cal_observe()

  elseif ARGS[1] == "2"
    scandisorder( ARGS[2], ARGS[3])

  elseif ARGS[1] == "3"
    truedisorder( ARGS[2], ARGS[3] )

  elseif ARGS[1] == "4"
    @time main(L=[3, 3], N=4, ex=20, int_ee=0.5, guess=false, sweepdim=200, sweepcnt=80)

  else
    println("not a valid operating mode")
  
  end 
end 

top()
#@time main(L=[3, 3], N=4, ex=20, guess=false)

