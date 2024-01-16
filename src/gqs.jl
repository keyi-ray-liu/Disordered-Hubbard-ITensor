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
        ψ_int = randomMPS( sites, state)

    elseif init_state == "emily_test_1"

        s1 = vcat([ "Occ" for _ in 1:N[2] ], [ "Emp" for _ in 1:N[1]])
        ψ1 = randomMPS( sites, s1)


        s2 = [ isodd(n) ? "Emp" : "Occ"  for n in 1:sum(N)]
        ψ2 = randomMPS( sites, s2)

        ψ_int = add( sqrt(0.9) * ψ1, sqrt(0.1) * ψ2)

    else
        error("init state string not defined")
    end 

    

    return ψ_int, sites
end


function GQS_measure(ψ)

    
    O = outer(ψ', ψ)
    s = siteinds(ψ)

    typestr = ["Emp", "Occ"]
    states = collect(Iterators.product(typestr, typestr))

    lefts = collect(1:div(length(ψ), 2))

    res = zeros( ComplexF64, (length(lefts), length(states)))


    for left in lefts
        for (col, state) in enumerate(states)

            basis = randomMPS( s[left:left+1], state)
            pO = partial_tr(O; leftend=left - 1, rightend=left + 2)
            res[left, col] = inner(basis', pO, basis)

        end 
    end

    return res

    #s = partial_contract(ψ, length(ψ), length(ψ) )

    # p = siteinds("S=1/2",12)

    # s = randomMPS(p)
    # s = partial_contract(s; leftend=100, rightend=100)
    # @show s

end 