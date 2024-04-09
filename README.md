# Disordered-Hubbard-ITensor
Code for doing 1D and small 2D simulations with ITensor DMRG package based on the ITensor algorithm.

## old
The top.jl contains the top level interface, as well as operation mode control.
The purpose of the individual files are self-explanatory.


## Current active code

The more Julian version of the code is currently in revamp-src. 
To run a specific program functionality, use

``` julia main.jl ARGUMENT ```

For example for NF calculations, 

``` julia main.jl NF_square ```

The parameters are loaded through JSON files, for example, to run a calculation with the following parameters:

``` {"U"=4.0, "L"=5, "Nup"=13, "Ndn"=11, "t"=0.001, "bias"=1.0} ```
Where "bias" is the value the center region will be biased with. 

You should write this JSON to a file "NFpara.json" under the main directory
