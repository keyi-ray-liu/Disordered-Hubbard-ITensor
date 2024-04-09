function SD_dynamics_QE(simu_para, sd_hop, additional_para)

    workdir = getworkdir()
    L = simu_para[:L]
    N = simu_para[:N]
  
    source_config = additional_para["source_config"]
    drain_config  = additional_para["drain_config"]
  
    Ntotal = N + count_ele(source_config) + count_ele(drain_config)
    Ltotal = get_systotal(para)
  
    if length(L) == 1
      sd_loc = [ [-1.0], [Ltotal]]
      source_site = 1
      drain_site = Ltotal
    else
      y = (minimum(L) - 1) / 2
      sd_loc = [ [y, -1.0 ], [y, maximum(L)]]
      source_site = div(L, 2) + 1
      drain_site = Ltotal - div(L, 2)
    end 
  
    sd_hop["sd_loc"] = sd_loc
    sd_hop["source_site"] = source_site
    sd_hop["drain_site"] = drain_site
  
    
  
    product_state = additional_para["product_state"]
    offset_output = "get_offset_bulk"
    initial_state_output = "initial_state"
    offset = sd_hop["source_offset"]
  
  
    τ = additional_para["τ"]
    start = additional_para["start"]
    fin = additional_para["fin"]
    occ_direct = additional_para["occ_direct"]
    
    if isnothing(offset) 
      # if no offset, we find the GS of the bulk
      if !isfile(workdir * offset_output * ".h5") 
        # if there no target file, we perform a single GS search, to find the correct offset for the SD
        println("No offset file, generating")
        paras = setpara(;simu_para..., N=Ntotal, ex=1, output = offset_output, s_len = 0, d_len = 0)
        main(paras; source_config = [], drain_config = [])
      else
        println("offset file found, loading")
      end 
  
      offset = readdlm(workdir * offset_output * "ex" )[end]
      sd_hop["source_offset"] = offset
      sd_hop["drain_offset"] = offset
    end 
  
    simu_para[:sd_hop] = sd_hop
  
    # we prepare the initial state
    if !product_state
      # if product state flag is off, we solve fo the GS of the bulk with initial configuration 
  
      if !isfile( workdir * initial_state_output * ".h5")
        # if no IS, generate one 
        paras = setpara(;simu_para..., ex=1, output = initial_state_output, sd_override=true)
        main(paras; source_config = source_config, drain_config = drain_config)
  
      end 
  
      wf = h5open( workdir * initial_state_output * ".h5", "r")
  
      if start == τ
        ψ = read(wf, "psi1", MPS)
      else 
        ψ = read(wf, "psi", MPS)
      end 
  
      sites = siteinds(ψ)
  
    else
      paras = setpara(;simu_para...)
      ψ, sites = TE_stateprep(paras; source_config = source_config, drain_config = drain_config )
    end 
  
    paras = setpara(;simu_para..., τ=τ,  output="TE")
    # further preparation of the initial state, if needed
    #ψ, sites = TE_stateprep(ψ, paras, sites)
    #ψ, sites = TE_stateprep(paras)
  
    println("length of ψ is:" , length(ψ))
    println("length of sites", length(sites))
  
    time_evolve(ψ, sites, paras, start, fin, occ_direct)
      
      
  
end 

"""
SD solver to solve the transport dyanmics. 

L-bulk sys-R

General workflow of the solver: simu parameters contains the core parameters, sd_hop has the SD parameters that OVERRIDES the sd_hop in simu dictionary. Additional parameters has parameters such as timesteps, etc. 

The get_mix_energy function controlled by the sd_hop flag, generates all the important structure of the mix basis, including ordering. All logic are wrapped in the function

if offset is nothing, the program runs with only the system to find the appropriate offset values. 

At the current stage, if the product_state flag is false, we then run the DMRG algorithm to find the ground state of our system with 0 offset.

When we have the initial state, either by loading existing WF or generating, we then run the time dynamics.


All the logic in SD should be decoupled from the bulk system, such that the transport dynamics can be separately specified by sd_hop

We default that offset must be provided, the offset calculations are removed for now
"""
function SD_dynamics_transport(simu_para, sd_hop, additional_para; kwargs...)

    workdir = getworkdir()
    L = simu_para[:L]
    geometry = simu_para[:geometry]


    source_config = additional_para["source_config"]
    drain_config  = additional_para["drain_config"]

    # N = simu_para[:N]
    # N_sys = N[2] + N[3] + 2 * N[4]
    # Ntotal = N_sys + count_ele(source_config) + count_ele(drain_config)

    Ltotal = get_systotal(L, geometry)
  
    if length(L) == 1
      sd_loc = [ [1.0], [Ltotal]]
      source_site = 1
      drain_site = Ltotal
    else
      y = (minimum(L) - 1) / 2
      sd_loc = [ [y, -1.0 ], [y, maximum(L)]]
      source_site = div( minimum(L) - 1, 2) + 1
      drain_site = Ltotal - div( minimum(L) - 1,  2)
    end 
  
    sd_hop["sd_loc"] = sd_loc
    sd_hop["source_site"] = source_site
    sd_hop["drain_site"] = drain_site
    
    bulk_bias = get(sd_hop, "bulk_bias", 0.0)
    init_bulk_bias = get(sd_hop, "init_bulk_bias", 0.0)

    if typeof(bulk_bias) == Float64 || typeof(bulk_bias) == Int
        bulk_bias = [bulk_bias for _ in 1:Ltotal]
    end 

    if typeof(init_bulk_bias) == Float64 || typeof(init_bulk_bias) == Int
      init_bulk_bias = [init_bulk_bias for _ in 1:Ltotal]
    end 

    if length(bulk_bias) != Ltotal || length(init_bulk_bias) != Ltotal
        error("bulk bias and sys size value mismatch")
    end
    
    sd_hop["bulk_bias"] = bulk_bias

    sdcouple = get(sd_hop, "sdcouple", 0)
    extends = []

    for extend in 1:sdcouple
      append!(extends, length(source_config) - extend + 1 )
      append!(extends, length(source_config) + extend + L)
    end 

    sd_hop["sdcouple"] = extends
    state_override = get(sd_hop, "state_override", [])
    mixbasis = sd_hop["mixbasis"]
    
    
    product_state = additional_para["product_state"]
    τ = additional_para["τ"]
    start = additional_para["start"]
    fin = additional_para["fin"]
    occ_direct = additional_para["occ_direct"]


    initial_state_output = "initial_state"
    
    # offset = sd_hop["source_offset"]
    # offset_output = "get_offset_bulk"
    # if isnothing(offset) 
    #   # if no offset, we find the GS of the bulk
    #   if !isfile(workdir * offset_output * ".h5") 
    #     # if there no target file, we perform a single GS search, to find the correct offset for the SD
    #     println("No offset file, generating")
    #     paras = setpara(;simu_para..., N=Ntotal, ex=1, output = offset_output, s_len =0, d_len = 0)
    #     main(paras; source_config = [], drain_config = [])
    #   else
    #     println("offset file found, loading")
    #   end 
  
    #   offset = readdlm(workdir * offset_output * "ex" )[end]
    #   sd_hop["source_offset"] = offset
    #   sd_hop["drain_offset"] = offset
    # end 
  
    simu_para[:sd_hop] = sd_hop

    # we will keep this basis even for the μ = 0 GS calculations for consistency
    if mixbasis
      mixbasis_energies, ks, LR, zeropoint = get_mix_energy(simu_para)
      addtags = "mix"

      println("zeropoint", zeropoint)
    else
      addtags = ""
      mixbasis_energies = ks = LR = []
    end 
  
    # we prepare the initial state
    if !product_state
      # if product state flag is off, we solve fo the GS of the bulk with initial configuration 
        

      # here we solve for the GS, where initially the chemical potential is turned off
      gs_para = deepcopy(simu_para)


        gs_para[ :sd_hop]["source_offset"] = gs_para[ :sd_hop]["drain_offset"] = 0
        gs_para[ :sd_hop]["bulk_bias"] = init_bulk_bias

      if mixbasis
        mixbasis_gs , _, _, _= get_mix_energy(gs_para; ks=ks, LR=LR)

      else
        mixbasis_gs = []
      end 

      if !isfile( workdir * initial_state_output * ".h5")
        # if no IS, generate one 
        paras = setpara(;gs_para..., ex=1, output = initial_state_output, sd_override=false)
        main(paras; mixbasis_energies = mixbasis_gs, ks= ks, LR =LR, addtags = addtags, source_config = source_config, drain_config = drain_config, state_override = state_override, kwargs...)
  
      end 
  
      wf = h5open( workdir * initial_state_output * ".h5", "r")
  
      if start == τ
        ψ = read(wf, "psi1", MPS)
      else 
        ψ = read(wf, "psi", MPS)
      end 
  
      sites = siteinds(ψ)
  
    else
      paras = setpara(;simu_para...)
      ψ, sites = TE_stateprep(paras; addtags=addtags, source_config = source_config, drain_config = drain_config, state_override = state_override)
    end 
  
    paras = setpara(;simu_para..., τ=τ,  output="TE")
    # further preparation of the initial state, if needed
    #ψ, sites = TE_stateprep(ψ, paras, sites)
    #ψ, sites = TE_stateprep(paras)
  
    println("length of ψ is:" , length(ψ))
    println("length of sites", length(sites))
    
    if !isempty( state_override)

      gs_n = expect(ψ, "N")
      mid = div(length(ψ), 2)

      println("temp occ check begin:")
      println( "sys:", gs_n[mid: mid+1])
      println( "sum of L:", sum(gs_n[begin:mid - 1]))
      println( "sum of R:", sum(gs_n[mid + 2:end]))

    end 

    if mixbasis

      #@show hastags(sites, "mix")
      writedlm(workdir * "mixbasis_energies", mixbasis_energies)
      writedlm( workdir * "mixbasis_gs", mixbasis_gs)
      writedlm( workdir * "ks", ks)
      writedlm( workdir * "LR", LR)
    end 

    time_evolve(ψ, sites, paras, start, fin, occ_direct; kwargs..., mixbasis_energies = mixbasis_energies, ks=ks, LR=LR)
      
    
  
end 


"""
Running TN simulations for localized transition. The sys N should be fixed at 2
"""

