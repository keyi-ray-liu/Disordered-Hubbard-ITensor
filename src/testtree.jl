using DataGraphs: edge_data, vertex_data
using Dictionaries: Dictionary
using Graphs: nv, vertices
using ITensors: ITensors
using ITensors.ITensorMPS: ITensorMPS
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
using ITensorNetworks.ModelHamiltonians: ModelHamiltonians
using KrylovKit: eigsolve
using NamedGraphs: named_comb_tree, rem_vertex!, add_vertex!, add_edge!, NamedGraph
using Observers: observer
using Test: @test, @test_broken, @testset
using ITensorUnicodePlots: @visualize


let 

  g = NamedGraph()
  L = 4
  ARM = 4

  for i in 1:ARM
    for j in 1:L
      add_vertex!(g, (i, j))
    end 
  end 

  add_vertex!(g, (5, 1))

  for i in 1:ARM

    add_edge!(g, ((5,1), (i, 1)))
    for j in 1:L - 1
      add_edge!(g, ((i, j), (i, j + 1)))
    end 
  end 

  

  @visualize g

end 



@testset "Tree DMRG for Fermions" for nsites in [2] #ToDo: change to [1,2] when random_ttn works with QNs
    auto_fermion_enabled = ITensors.using_auto_fermion()
    use_qns = true
    cutoff = 1e-12
    nsweeps = 10
    maxdim = [10, 20, 40, 100]
  
    # setup model
    tooth_lengths = fill(2, 3)
    c = named_comb_tree(tooth_lengths)

    s = siteinds("Electron", c; conserve_qns=use_qns)
    U = 2.0
    t = 1.3
    tp = 0.6
    os = ModelHamiltonians.hubbard(c; U, t, tp)
  
    # for conversion to ITensors.MPO
    linear_order = [4, 1, 2, 5, 3, 6]
    vmap = Dictionary(vertices(s)[linear_order], 1:length(linear_order))
    sline = only.(collect(vertex_data(s)))[linear_order]
  
    # get MPS / MPO with JW string result
    ITensors.disable_auto_fermion()
    Hline = ITensorMPS.MPO(relabel_sites(os, vmap), sline)
    psiline = ITensorMPS.randomMPS(sline, i -> isodd(i) ? "Up" : "Dn"; linkdims=20)
    e_jw, psi_jw = dmrg(Hline, psiline; nsweeps, maxdim, cutoff, outputlevel=0)
    ITensors.enable_auto_fermion()
  
    # now get auto-fermion results 
    H = ttn(os, s)
    # make init_state
    d = Dict()
    for (i, v) in enumerate(vertices(s))
      d[v] = isodd(i) ? "Up" : "Dn"
    end
    states = v -> d[v]
    psi = ttn(s, states)
    psi = dmrg(
      H, psi; nsweeps, maxdim, cutoff, nsites, updater_kwargs=(; krylovdim=3, maxiter=1)
    )
  
    # Compare to `ITensors.MPO` version of `dmrg`
    Hline = ITensorMPS.MPO(relabel_sites(os, vmap), sline)
    psiline = ITensorMPS.randomMPS(sline, i -> isodd(i) ? "Up" : "Dn"; linkdims=20)
    e2, psi2 = dmrg(Hline, psiline; nsweeps, maxdim, cutoff, outputlevel=0)
  
    @test inner(psi', H, psi) ≈ inner(psi2', Hline, psi2) atol = 1e-5
    @test e2 ≈ e_jw atol = 1e-5
    @test inner(psi2', Hline, psi2) ≈ e_jw atol = 1e-5
  
    if !auto_fermion_enabled
      ITensors.disable_auto_fermion()
    end
  end