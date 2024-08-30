
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


const TEMP_tag = "temp_temp_"

@show Regex(TEMP_tag * get_static_str("biasedchain") * "*.h5")
occursin(Regex(TEMP_tag * get_static_str("biasedchain") * ".*.h5"), "temp_temp_plasmon1111111.h5")