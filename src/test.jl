

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
