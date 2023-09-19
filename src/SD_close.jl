function SD_dynamics(simu_para, sd_hop, additional_para)

    workdir = getworkdir()
    L = simu_para[:L]
    N = simu_para[:N]
  
    source_config = simu_para[:source_config]
    drain_config  = simu_para[:drain_config]
  
    Ntotal = N + count_ele(source_config) + count_ele(drain_config)
    Ltotal = prod(L)
  
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
        paras = setpara(;simu_para..., N=Ntotal, ex=1, output = offset_output, source_config = [], drain_config = [])
        main(paras;)
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
        main(paras;)
  
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
      ψ, sites = TE_stateprep(paras)
    end 
  
    paras = setpara(;simu_para..., τ=τ,  output="TE")
    # further preparation of the initial state, if needed
    #ψ, sites = TE_stateprep(ψ, paras, sites)
    #ψ, sites = TE_stateprep(paras)
  
    println("length of ψ is:" , length(ψ))
    println("length of sites", length(sites))
  
    time_evolve(ψ, sites, paras, start, fin, occ_direct)
      
      
  
end 


function SD_dynamics_transport(simu_para, sd_hop, additional_para)

    workdir = getworkdir()
    L = simu_para[:L]
    N = simu_para[:N]
  
    source_config = simu_para[:source_config]
    drain_config  = simu_para[:drain_config]
  
    Ntotal = N + count_ele(source_config) + count_ele(drain_config)
    Ltotal = prod(L)
  
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
        paras = setpara(;simu_para..., N=Ntotal, ex=1, output = offset_output, source_config = [], drain_config = [])
        main(paras;)
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
        

      # here we solve for the GS, where initially the chemical potential is turned off
      gs_para = deepcopy(simu_para)
      gs_para[ :sd_hop]["source_offset"] = gs_para[ :sd_hop]["drain_offset"] = 0

      if !isfile( workdir * initial_state_output * ".h5")
        # if no IS, generate one 
        paras = setpara(;gs_para..., ex=1, output = initial_state_output, sd_override=false)
        main(paras;)
  
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
      ψ, sites = TE_stateprep(paras)
    end 
  
    paras = setpara(;simu_para..., τ=τ,  output="TE")
    # further preparation of the initial state, if needed
    #ψ, sites = TE_stateprep(ψ, paras, sites)
    #ψ, sites = TE_stateprep(paras)
  
    println("length of ψ is:" , length(ψ))
    println("length of sites", length(sites))
  
    time_evolve(ψ, sites, paras, start, fin, occ_direct)
      
      
  
end 

