function GQS_dynamic(gqs_para, additional_para)

    #TE para
    init_state = additional_para["init_state"]

    τ = additional_para["τ"]
    start = additional_para["start"]
    fin = additional_para["fin"]
    occ_direct = additional_para["occ_direct"]
    
    para = setpara(;gqs_para...)
    ψ_int, sites = GQS_init(para, init_state)

    para = setpara(;gqs_para..., τ=τ, output = "TE")
    time_evolve(ψ_int, sites, para, start, fin, occ_direct)
end


function GQS_init(para, init_state)

    sites = init_site(para)

    N = para["N"]

    if init_state == "all_left_spinless"
        state = vcat([ "Occ" for _ in 1:N[2] ], [ "Emp" for _ in 1:N[1]])

    else
        error("init state string not defined")
    end 

    ψ_int = randomMPS( sites, state)

    return ψ_int, sites
end
