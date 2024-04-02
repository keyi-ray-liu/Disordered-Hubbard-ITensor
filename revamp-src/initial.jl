gen_state_str(sys::systems) = shuffle([ get_type_dict(type(sys))[i] for i= 1:4 for _ in 1:N(sys)[i] ])

gen_state_str(sys::QE_two) = vcat(QE_str(sys)[1], gen_state_str(sys.chain_only), QE_str(sys)[2])


gen_state_str(sys::DPT) = vcat([ get_type_dict(type(sys))[i] for i=1:4 for _ in 1:N(sys)[i]],  [ get_type_dict(type(sys))[i] for i=1:4 for _ in 1:N(sys)[i]], ["Occ", "Emp"])

gen_state_str(sys::LSR_SIAM) = vcat([ get_type_dict(type(sys))[i] for i=1:4 for _ in 1:N(sys)[i]], ["Emp"], [ get_type_dict(type(sys))[i] for i=1:4 for _ in 1:N(sys)[i]])


function QE_str(sys::QE_two)

    ex = ["Occ", "Emp"]
    gs = ["Emp", "Occ"]

    if init(sys) == "Left"
        return ex, gs

    elseif init(sys) == "Right"
        return gs, ex

    elseif init(sys) == "Both"
        return ex, ex

    elseif init(sys) == "None"
        return gs, gs

    else
        error("Unrecognized QE init")
    end 

end 

function gen_state(sys::systems; QN=true, kwargs...)


    sites = siteinds(type(sys), get_systotal(sys); conserve_qns=QN)

    #@show length(gen_state_str(sys)), length(get_systotal(sys))
    @show gen_state_str(sys)
    ψ = randomMPS(sites, gen_state_str(sys) )

    return ψ

end 



function gen_hamiltonian(sys::systems)

    res = OpSum()

    res = add_hop!(sys, res)
    res = add_DensityDensity!(sys, res)
    res = add_onsite!(sys, res)
    res = add_qe!(sys, res)
    res = add_HubbardRepulsion!(sys, res)
    res = add_specific_int!(sys, res)
    
    return res
end 


# All the system specific logic is wrapped with multiple-dispatch neighbor functions, we no longer perform these logics within the code





