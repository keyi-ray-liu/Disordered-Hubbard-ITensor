

function runtest()
    N = 200

    s = siteinds("Fermion", N;
    conserve_qns=true
    )
    
    state = [ isodd(n) ? "Occ" : "Emp" for n in 1:N]
    M = randomMPS(s, state )

    
    ψ = load_ψ( "temp_cur_ex1"; tag="psi")
    L = length(ψ)

    M[1:L] = ψ

    typeof(M)
    @show M

    @show expect(M, "N")
    return nothing

end 

function corr_test()

    qn = true

    L = 80
    N = 20
    
    t = -1.0
    s = siteinds("Electron", L;
    conserve_qns=qn
    )
    
    @show leftinds = shuffle([isodd(n) ? 1 : 0 for n in 1:L])
    @show latticeleftind = positiveind(leftinds)

    # Free fermion hopping Hamiltonian
    #h = Hermitian(diagm(1 => fill(-t, N-1), -1 => fill(-t, N-1)))

    h = zeros(L, L)

    for i in eachindex(latticeleftind[1:end - 1])
        
        h[ latticeleftind[i], latticeleftind[i + 1]] = t
        h[ latticeleftind[i + 1], latticeleftind[i]] = t
        h[ latticeleftind[i], latticeleftind[i]] = -10

    end 

    _, u = eigen(h)

    # Get the Slater determinant
    Φ = u[:, 1:N]
    ψ0 = slater_determinant_to_mps(s, Φ, Φ; maxblocksize = 4)

    
    @show expect(ψ0, "Ntot")


end 


function init_test()

    N = 40
    L = 20
    
    @show state = [  "Emp" for _ in 1:N]

    s = siteinds("Electron", N; conserve_qns=true)
    
    #ref = [ n <= L ? "Emp" : "Occ" for n in 1:N]

    M = randomMPS(s, state)

    for j in L:1

        # operator = ITensor(1.0)
        # for k in 1:N

        #     U = Ujk(N, j, k)
        #     operator *= U * op("Cdag", s[k])
            
        # end 
        operator = op("Cdagup", s[j])
        M = apply(operator, M)
        operator = op("Cdagdn", s[j])
        M = apply(operator, M)

        @show j
        operator = op("Cdagdn", s[j])
        M = apply(operator, M)
    end 

    # Ms = [ sqrt(A[i]) * randomMPS(s, [ n == i ? "Occ" : "Emp" for n in 1:N] ) for i in 1:40]
    
    # M = add( Ms..., maxdim=10)

    # M = add(M, reference)

    @show expect(M, "Ntot")
    return nothing


end 

#corr_test()

# @show Regex(TEMP_tag * get_static_str("biasedchain") * "*.h5")
# occursin(Regex(TEMP_tag * get_static_str("biasedchain") * ".*.h5"), "temp_temp_plasmon1111111.h5")


function corr_test2()

    L = 4
    s = siteinds("Electron", L; conserve_qns = true)

    str = [ isodd(i) ? "Up" : "Dn" for i in 1:L]

    M = randomMPS(s, str)

    @show (correlation_matrix(M, "Sz", "Sz"))


end 

function typetest()

    sys = Chain()
    @show  dis(0, 1, sys)
    res = []
    @show add_specific_int!(sys, res)

end 

function argtest(; kwargs...)

    L = get(kwargs, :L, error("No L!"))


end 






function lapacktest()

    N = 128
    bias = -1e3
    sys = 9
    total = N * 2 + sys

    s = siteinds("Electron", total; conserve_qns = true)
    state = [ i<= N ? "UpDn" : "Emp" for i in 1:total]

    M = randomMPS(s, state)

    H = OpSum()

    for i in 1:N - 1
        H += -1.0, "Cdagup", i, "Cup", i + 1
        H += -1.0, "Cdagup", i + 1, "Cup", i
        H += -1.0, "Cdagdn", i, "Cdn", i + 1
        H += -1.0, "Cdagdn", i + 1, "Cdn", i
    end 

    for i in N + sys + 1: total - 1
        H += -1.0, "Cdagup", i, "Cup", i + 1
        H += -1.0, "Cdagup", i + 1, "Cup", i
        H += -1.0, "Cdagdn", i, "Cdn", i + 1
        H += -1.0, "Cdagdn", i + 1, "Cdn", i
    end 


    sc = -0.25
    H += sc, "Cdagup", N, "Cup", N + 1
    H += sc, "Cdagup", N + 1, "Cup", N
    H += sc, "Cdagdn", N, "Cdn", N + 1
    H += sc, "Cdagdn", N + 1, "Cdn", N

    d = -0.25
    H += d, "Cdagup", N + sys, "Cup", N + sys + 1
    H += d, "Cdagup", N + sys + 1, "Cup", N + sys
    H += d, "Cdagdn", N + sys, "Cdn", N + sys + 1
    H += d, "Cdagdn", N + sys + 1, "Cdn", N + sys

    for i in N + 1:N+ sys

        r = (i - N - 1) ÷ 3 + 1
        c = (i - N - 1) % 3 + 1
        tos = []

        if r < 3
            append!(tos, [i + 3])
        end 

        if c < 3
            append!(tos, [i + 1])
        end 

        @show i, tos

        for to in tos
            H += -1, "Cdagup", i, "Cup", to
            H += -1, "Cdagup", to, "Cup", i
            H += -1, "Cdagdn", i, "Cdn", to
            H += -1, "Cdagdn", to, "Cdn", i
        end 

        sumval = 0
        for j in (i + 1) : (N + sys)

            rj = (j - N - 1) ÷ 3 + 1
            cj = (j - N - 1) % 3 + 1

            H += 1/( sqrt( (rj - r) ^2 + (cj - c)^2) + 0.5), "Ntot", i, "Ntot", j

            sumval += -1/( sqrt( (rj - r) ^2 + (cj - c)^2) + 0.5)
        end 


        H += 4.0, "Nupdn", i
        H += sumval, "Ntot", i
    end 


    for i in 1:N
        H += bias, "Ntot", i 
    end 

    for i in N + sys + 1: total
        H += - bias, "Ntot", i 
    end 


    H = MPO(H, s)

    sweep = Sweeps(20)
    setnoise!(sweep, 1E-6)
    setmaxdim!(sweep, 128)
    dmrg(H, M, sweep)
end 


function toytwolevel()

    L = 32
    total = L * 2 + 2
    μC = -1.23
    vs = 0.25

    dd = 1.3
    s = siteinds("Fermion", total; conserve_qns = true)

    ifshuffle = false
    f(x) = ifshuffle ? shuffle(x) : x
    state = f([ isodd(i) ? "Occ" : "Emp" for i in 1:total])

    M = randomMPS(s, state)

    H = OpSum()

    for i in L - 1: L + 2
        H += dd, "N", i
    end 

    for i in 1:total - 2 - 1
        H += 1, "Cdag", i, "C", i + 1
        H += 1, "Cdag", i + 1, "C", i
    end 

    H += μC, "N", total - 1
    H += vs, "Cdag", total - 1, "C", total
    H += vs, "Cdag", total, "C", total - 1


    saveham("toy", H)
    H = MPO(H, s)

    

    sweep = Sweeps(20)
    setnoise!(sweep, 1E-6)
    setmaxdim!(sweep, 64)
    w, ψ =  dmrg(H, M, sweep)

    @show correlation_matrix(ψ, "Cdag", "C", sites= total - 1 : total)


end 



function test_rotation_gate()

    s = siteinds("Fermion", 2; conserve_qns = true)
    state1 = ["Emp", "Emp"]
    #state2 = ["Emp", "Occ"]



    #a += sqrt(0.7), "Cdag", 1
    #a += sqrt(0.3), "Cdag", 2

    #a = sqrt(0.7) * Op("cdag", s, 1)
    #b = sqrt(0.3) * Op("cdag", s, 2)


    # a = sqrt(0.7) * op("Cdag", s[1]) 
    # @show a


    M = zeros(ComplexF64, 4, 4)

    M[2, 1] = sqrt(0.6)
    M[3, 1] = sqrt(0.4) * exp( 1im * 0.5)
    #sqrt(0.7)
    #M[2, 2] = sqrt(0.3)

    a = op(M, s[1], s[2])

    ψ = randomMPS(s, state1)
    ψ = apply(a, ψ)
    #ψ2 = randomMPS(s, state2)

    #ψ :: MPS = add(sqrt(0.7) * ψ1, sqrt(0.3) * ψ2)

    # @show typeof(ψ), fieldnames(typeof(ψ))
    # @show typeof(ψ.data), fieldnames(typeof(ψ.data))
    # @show typeof(ψ.data[1]), fieldnames(typeof(ψ.data[1]))
    # @show fieldnames(typeof(ψ.data[1].tensor.storage))
    #@show typeof(ψ.data[1].tensor.inds), fieldnames(typeof(ψ.data[1].tensor.inds))




    #ψ.data[1].tensor.storage.data = [sqrt(0.1), sqrt(0.9)]
    @show expect(ψ, "N")
    @show correlation_matrix(ψ, "Cdag", "C")

end 


function test_hop_gate()

    s = siteinds("Fermion", 2; conserve_qns = true)
    state1 = ["Occ", "Emp"]
    #state2 = ["Emp", "Occ"]
    
    a = OpSum()
    a += sqrt(0.7), "Cdag", 1, "C", 2
    a += sqrt(0.3), "Cdag", 2, "C", 1

    a = MPO(a, s)
    ψ = randomMPS(s, state1)
    ψ = apply(a, ψ)

    @show expect(ψ, "N")
    @show correlation_matrix(ψ, "Cdag", "C")

end 



# function test_DPT()


    

#     dpt_in = Dict()


#     L = get(dpt_in, "L", 34)
#     R = get(dpt_in, "R", 34)
#     tswitch = Float64(get(dpt_in, "tswitch", 0.0))
#     t_fin = Float64(get(dpt_in, "tfin", 30.0))
#     τ = get(dpt_in, "timestep", 0.25)
#     TEdim = get(dpt_in, "TEdim", 64)
#     biasLR = get(dpt_in, "biasLR", 0.0)
#     mixed = get(dpt_in, "mixed", true)
#     ordering = get(dpt_in, "ordering", "SORTED")
#     avg = get(dpt_in, "avg", false)
#     sweepcnt = get(dpt_in, "sweepcnt", 30)
#     vs = get(dpt_in, "vs", 0.25)

#     κs = [sqrt(0.5), sqrt(0.8), sqrt(0.95)]
#     θs = [0, 1/3 * π, 2/3 * π, π]
#     γs = [0.1, 0.5, 1.0]
#     Us = [0.1, 3.25, 5.5]

#     Threads.@threads for (κ, θ, γ, U) in collect(Base.product(κs, θs, γs, Us ))

#         ψ1 = run_DPT_many_body(U, L, R,  0.0; tswitch = 0.0, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  ddposition="R",  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "manuallower", vs = vs, ordering = ordering, workflag = "lok$(κ)t$(θ)g$(γ)U$(U)")

#         sites = siteinds(ψ1)

#         ψ2 = run_DPT_many_body(U, L, R,  0.0; tswitch = 0.0, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  ddposition="R",  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "manualupper", sites = sites, vs = vs, ordering = ordering, workflag = "upk$(κ)t$(θ)g$(γ)U$(U)")

#         ψ = add( κ * ψ1, γ * exp(1im* θ) * sqrt(1 - κ^2) * ψ2)


#         expval = expect(ψ, "N")
#         corr = correlation_matrix(ψ, "Cdag", "C")

#         @show expval[end - 1: end]
#         @show corr[end - 1: end, end - 1: end]

#         process = Supplywf(ψ)

#         run_DPT_many_body(U, L, R,  t_fin; tswitch = tswitch, bias_L = biasLR/2, bias_R  = - biasLR/2, τ=τ, mixed=mixed,  workflag = "_k$(κ)_t$(θ)_g$(γ)_U$(U)", ddposition="R",  avg=avg,   TEdim = TEdim, sweepcnt = sweepcnt, mode = "override", vs = vs, process = process)

#     end 

# end 



function test_corr_MPO()


    N = 4
    s = siteinds("Electron", N; conserve_qns = true )
    state = [ isodd(i) ? "Emp" : "Up" for i in 1:N]
    psi1 = randomMPS(s, state)

    state2 = [ iseven(i) ? "Emp" : "Up" for i in 1:N]
    psi2 = randomMPS(s, state2)

    psi = add( sqrt(0.7) * psi1, sqrt(0.3) * psi2)
    @show expect(psi, "Sz")

    a = OpSum()
    for i in 1:N
        for j in 1:N
            a += 1.0, "Sz", i, "Sz", j
        end 
    end 

    A = MPO(a, s)

    @show inner(psi', A, psi)
    @show sum(correlation_matrix(psi, "Sz", "Sz"))

end 