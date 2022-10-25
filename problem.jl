using ITensors
using ITensorTDVP: linsolve

function test()
    L = 12
    N = 6
    t = 1.0

    λ = - 10.0
    sites = siteinds("Fermion", L, conserve_qns=true)
    ampo = OpSum()

    for i=1:L-1
        ampo += -t , "C",i,"Cdag",i+1
        ampo += -t , "C",i+1,"Cdag",i
    end 

    H = MPO(ampo,sites)
    state = append!([ "Occ" for n=1:N] , ["Emp" for n=1:L -N])
    ϕ = randomMPS(sites,state)

    # outer iteration, solving excited state eigen energy immediately above λ
    for ex=1:10

        # power iteration, repeatedly solving (H - λ) ψ = ϕ, ψ → ϕ
        for cnt = 1:100

            # generate an initial guess ψ0
            #state = append!([ "Occ" for n=1:N] , ["Emp" for n=1:L -N])
            #ψ0 = randomMPS(sites,state)
            ψ0 = ϕ
            ψ = linsolve(H, ϕ, ψ0, -λ, 1.0)
            ϕ = ψ

            energy = inner(ϕ', H, ϕ) / inner( ϕ', ϕ)
            println("energy = ", energy)

        end 

        λ = energy
        println("λ =  ", λ)
    end 

end 

test()
