add_penalty!(sys :: Systems, res :: OpSum ) = res

add_penalty!(sys :: DPT_mixed, res :: OpSum) = add_penalty!(sys.dpt, res)

function add_penalty!( sys:: DPT, res :: OpSum)

    # we add a penalty of the form V(n - n0)^2
    if !isnothing(sys.n1penalty)

        res += 1000, "N", dd_lower(sys), "N", dd_lower(sys)
        res -= 1000 * 2 * sys.n1penalty, "N", dd_lower(sys)
    end 

    return res
end 