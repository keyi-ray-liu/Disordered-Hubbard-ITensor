

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
    N = 40
    L = 20
    

    s = siteinds("Fermion", N;
    conserve_qns=qn
    )

    if !qn

        @show A = rand(N)
        M = randomMPS(s)
        
        for i in 1:N
            a = A[i]
            M[i] = ITensor([sqrt(1 -a), sqrt(a)], s[i])
        end
    else

        @show state = [ n <= L ? "Occ" : "Emp" for n in 1:N]

        M = randomMPS(s, state)

        for i in 1:N
            a = 0.5 #A[i]
            M[i] = ITensor([sqrt(1 -a), sqrt(a)], s[i])
        end

    end

    # Ms = [ sqrt(A[i]) * randomMPS(s, [ n == i ? "Occ" : "Emp" for n in 1:N] ) for i in 1:40]
    
    # M = add( Ms..., maxdim=10)

    @show expect(M, "N")
    return nothing


end 

#corr_test()

# @show Regex(TEMP_tag * get_static_str("biasedchain") * "*.h5")
# occursin(Regex(TEMP_tag * get_static_str("biasedchain") * ".*.h5"), "temp_temp_plasmon1111111.h5")