
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

    para = Dict(
        :L => 100,
        :N => 25,
        :sweepdim => 400,
        :sweepcnt => 30,
        :QE => 2,
        :QEen => 0.1030502,
        :dynamode => "left",
        :TEmethod => "TDVP",
        :TEcutoff => 1E-9,
        :type => "Fermion",
        :krylovdim => 8,
        :QN=>true,
        :range=>1000,
        :range_qe => 1000
    )

    additional_para = Dict(
        "τ" => 0.5,
        "start" => 0.5,
        "fin" => 100.0,
        "product_state" => false,
        "occ_direct" => false,
        "QEmul" => 1.0
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
    additional_paras= JSON.parsefile( cur_dir * "/NFpara.json")

    paras = Dict{Any, Any}(
        :t => get(additional_paras, "t", -1.0),
        :N => get(additional_paras, "N", 0) ,
        :geometry => get(additional_paras, "geometry", "linear"),
        :L => get(additional_paras, "L", 0),
        :sweepdim => get(additional_paras, "sweepdim", 500),
        :sweepcnt => get(additional_paras, "sweepcnt", 100),
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
    transport_para = JSON.parsefile( cur_dir * "/transportpara.json")
    
    para = Dict{Any, Any}(
        :L => get(transport_para, "L", 1),
        :N => get(transport_para, "N", 0),
        :sweepdim => get(transport_para, "sweepdim", 128),
        :sweepcnt => get(transport_para, "sweepcnt", 150),
        :QEen => 0.0,
        :TEmethod => "TDVP",
        :TEcutoff => get(transport_para, "TEcutoff", 1E-9),
        :TEdim => get(transport_para, "TEdim", 128),
        :type => get(transport_para, "type", "Electron"),
        :krylovdim => get(transport_para, "krylovdim", 128),
        :QN=>true,
        :CN=> get(transport_para, "CN", 0),
        :int_ee => get(transport_para, "int_ee", 0.0),
        :int_ne => get(transport_para, "int_ne", 0.0),
        :U => get(transport_para, "U", 0.0),
        :t => get(transport_para, "t", 1.0),
    )
    
    sd_hop = JSON.parsefile( cur_dir * "/sdpara.json")

    additional_para = Dict(
        "τ" => 1.0,
        "start" => 1.0,
        "fin" => 100,
        "product_state" => false,
        "occ_direct" => false
    )

    occ1 = 1 + (para[:type] == "Fermion" ? 0 : 1)
    occ2 = occ1 + 1

    NR = get(sd_hop, "NR", 128)
    para[:source_config] = [ x in StatsBase.sample(1:NR, div(NR, 2), replace = false) ? occ1 : occ2 for x in 1:NR]
    para[:drain_config] = [ x in StatsBase.sample(1:NR, div(NR, 2), replace = false) ? occ1 : occ2 for x in 1:NR]

    #para[:source_config] = [ 2 for _ in 1:NR]
    #para[:drain_config] = [ 3 for _ in 1:NR]
    SD_dynamics_transport(para, sd_hop, additional_para)


end 


function time_obs_wrapper(num, obs)


    num = parse(Int, num)
    if num == 1
        method = "TEBD"

    elseif num == 2
        method = "TDVP"
    end 

    para = Dict{Any, Any}(
        "obs" => obs,
        "method" => method,
        "system" => "QE",
        "type" => "Fermion"
    )



    if obs == "tcd"
        para["λ_ee"] = 2.0
        para["ζ"] = 0.5
        para["section"] = 10
    end 

    if obs == "current"
        cur_dir = pwd()
        sd_hop = JSON.parsefile( cur_dir * "/sdpara.json")
        
        para["N"] = sd_hop["NR"]
        para["v"] = sd_hop["to_chain_hop"]
    end 

    time_obs(para)
end


