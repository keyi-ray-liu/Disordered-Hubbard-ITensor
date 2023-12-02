
function GS_wrapper()

    para = Dict(
        :L => 120,
        :N => 6,
        :ex => 2,
        :sweepdim => 300,
        :sweepcnt => 50,
        :krylovdim => 8,
        :range => 10000
    )

    GSGap(para)
end 


function QE_wrapper(num, energy)

    num = parse(Int, num)
    energy = parse(Float64, energy)

    para = Dict(
        :L => 100,
        :sweepdim => 400,
        :sweepcnt => 30,
        :QE => num,
        :QEen => energy,
        :ex => 20,
        :krylovdim => 8,
        :range=>1000,
        :range_qe => 1000
    )

    QE(para)
end 


function QEdyna_wrapper()
    cur_dir = pwd()


    #qedyna_string = read( cur_dir * "/QEdynapara.json", String)
    qedyna_paras = load_JSON(cur_dir * "/QEdynapara.json")

    para = Dict(
        :L => get(qedyna_paras, "L", 100),
        :N => get(qedyna_paras, "N", 50),
        :sweepdim => get(qedyna_paras, "sweepdim", 400),
        :sweepcnt => get(qedyna_paras, "sweepcnt", 30),
        :QE => get(qedyna_paras, "QE", 2),
        :QEen => get(qedyna_paras, "QEen", nothing),
        :dynamode => get(qedyna_paras, "dynamode", "left"),
        :TEmethod => get(qedyna_paras, "TEmethod", "TDVP"),
        :TEdim => get(qedyna_paras, "TEdim", 300),
        :TEcutoff => get(qedyna_paras, "TEcutoff", 1E-9),
        :type => get(qedyna_paras, "type", "Fermion"),
        :krylovdim => get(qedyna_paras, "krylovdim", 8),
        :QN=> get(qedyna_paras, "QN", true),
        :range=>get(qedyna_paras, "range", 1000),
        :range_qe => get(qedyna_paras, "range_qe", 1000),
        :int_ne => get(qedyna_paras, "int_ne", 2.0),
        :int_ee => get(qedyna_paras, "int_ee", 2.0),
        :exch => get(qedyna_paras, "exch", 0.2)
    )

    
    add_para = load_JSON(cur_dir * "/addpara.json")

    additional_para = Dict(
        "τ" => get(add_para, "t", 0.5),
        "start" => get(add_para, "start", 0.5),
        "fin" => get(add_para, "fin", 100.0),
        "product_state" => get(add_para, "product_state", false),
        "occ_direct" => get(add_para, "occ_direct", false),
        "QEmul" => get(add_para, "QEmul", 1.0)
    )

    QE_dynamic(para, additional_para)
end 


function eigensolver_wrapper()

    GS_para = Dict{Any, Any}(
        :L => 120,
        :ex => 2,
        :sweepdim => 300,
        :sweepcnt => 25,
        :krylovdim => 8,
        :QN => true,
        :range => 10000
      )

    QE_internal_para = Dict(
        "QE_mul" => 1,
        "pl_level" => 2
    )

    QE_para = deepcopy(GS_para)

    QE_para[:ex] = 4
    QE_para[:QE] = 2
    QE_para[:screening_qe] = 0.0

    dyna_para = deepcopy(QE_para)
    dyna_para[:dynamode]= "left"
    dyna_para[:TEmethod] = "eigen-occ"

    additional_para = Dict(
        "τ" => 1.0,
        "start" => 1.0,
        "fin" => 1000.0,
        "product_state" => false,
        "occ_direct" => false
    )

    println(GS_para)
    println(QE_para)
    println(dyna_para)
    
    eigensolver(GS_para, QE_internal_para, QE_para, dyna_para, additional_para)
end 


function corr_wrapper()

    para = Dict{Any, Any}(
        "op1" => "Cdag",
        "op2" => "C",
        "tag" => "CC",
        "wftag" => "tTDVP", 
    )

    time_corr_plot(para)
end


#function NF_wrapper(t, spdim, dim, Nup, Ndn, geometry)
function NF_wrapper()


    cur_dir = pwd()
    additional_paras= load_JSON( cur_dir * "/NFpara.json")

    paras = Dict{Any, Any}(
        :t => get(additional_paras, "t", -1.0),
        :N => get(additional_paras, "N", 0) ,
        :geometry => get(additional_paras, "geometry", "linear"),
        :L => get(additional_paras, "L", 0),
        :sweepdim => get(additional_paras, "sweepdim", 500),
        :sweepcnt => get(additional_paras, "sweepcnt", 150),
        :ex => 1,
        :krylovdim => get(additional_paras, "krylovdim", 10),
        :U => 4.0,
        :type => "Electron",
        :int_ee => 0.0,
        :int_ne => 0.0,
        :cutoff => get(additional_paras, "cutoff", 1E-16)
    )


    NF(paras)

end 


function iter_sd_wrapper()

    for mul in 1:20

        offset = mul * 0.2
        source_drain_wrapper(;offset=offset)

    end 

end 


function source_drain_wrapper(; offset=nothing)

    para = Dict{Any, Any}(
        :L => 12,
        :N => 1,
        :sweepdim => 100,
        :sweepcnt => 20,
        :QEen => 0.0,
        :dynamode => "left",
        :TEmethod => "TDVP",
        :TEcutoff => 1E-30,
        :type => "Fermion",
        :krylovdim => 8,
        :QN=>true,
        :CN=>1
    )
    
    sd_hop = Dict{Any, Any}(
        "to_chain_hop" => 0.001,
        "internal_hop" => 1.0,
        "source_offset" => offset,
        "drain_offset" => offset,
        "bulk_bias" => 0.0
    )

    additional_para = Dict(
        "τ" => 0.1,
        "start" => 0.1,
        "fin" => 10,
        "product_state" => true,
        "occ_direct" => true
    )

    para[:source_config] = [2]
    para[:drain_config] = [1]
    SD_dynamics(para, sd_hop, additional_para)


end 

function transport_wrapper()

    cur_dir = pwd()
    transport_para = load_JSON(cur_dir * "/transportpara.json")
    
    para = Dict{Any, Any}(
        :L => get(transport_para, "L", 1),
        :N => get(transport_para, "N", 0),
        :sweepdim => get(transport_para, "sweepdim", 128),
        :sweepcnt => get(transport_para, "sweepcnt", 150),
        :QEen => 0.0,
        :TEmethod => "TDVP",
        :cutoff => get(transport_para, "cutoff", 1E-10),
        :TEcutoff => get(transport_para, "TEcutoff", 1E-10),
        :TEdim => get(transport_para, "TEdim", 128),
        :type => get(transport_para, "type", "Fermion"),
        :krylovdim => get(transport_para, "krylovdim", 128),
        :QN=>true,
        :CN=> get(transport_para, "CN", 0),
        :int_ee => get(transport_para, "int_ee", 0.0),
        :int_ne => get(transport_para, "int_ne", 0.0),
        :U => get(transport_para, "U", 0.0),
        :t => get(transport_para, "t", 1.0),
    )
    
    sd_hop = load_JSON(cur_dir * "/sdpara.json")
    add_para = load_JSON(cur_dir * "/addpara.json")

    additional_para = Dict{String, Any}(
        "τ" => get(add_para, "t", 0.5),
        "start" => get(add_para, "start", 0.5),
        "fin" => get(add_para, "fin", 100.0),
        "product_state" => get(add_para, "product_state", false),
        "occ_direct" => get(add_para, "occ_direct", false)
    )


    source_config = sd_gen(sd_hop, which="source", type=para[:type])
    drain_config = sd_gen(sd_hop, which ="drain", type=para[:type])

    para[:s_len] = length(source_config)
    para[:d_len] = length(drain_config)

    additional_para["source_config"] = source_config
    additional_para["drain_config"] = drain_config

    #para[:source_config] = [ 2 for _ in 1:NR]
    #para[:drain_config] = [ 3 for _ in 1:NR]
    SD_dynamics_transport(para, sd_hop, additional_para)


end 


function time_obs_wrapper( obs)

    cur_dir = pwd()
    obs_para = load_JSON(cur_dir * "/obspara.json")

    para = Dict{Any, Any}(
        "obs" => obs,
        "method" => get(obs_para, "method", "TDVP"),
        "system" => get(obs_para, "system", "QE"),
        "type" => get(obs_para, "type", "Fermion"),
    )

    if obs == "tcd"
        para["λ_ee"] = 2.0
        para["ζ"] = 0.5
        para["section"] = 10
    end 

    if obs == "current"
        

        sd_hop = load_JSON(cur_dir * "/sdpara.json")
        para["N"] = sd_hop["NR"]
        para["v"] = sd_hop["to_chain_hop"]

    end 

    time_obs(para)
end


