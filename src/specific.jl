

Uk(j::Int, sys::DPT_mixed) = Uk(j, ks(sys), LR(sys))

"""generate the  transformation LR vectors for the reservoir(s) SPATIAL site j, with given ordering of the momentum space, with LR cross terms being zero"""
function Uk(j::Int, ks, LR)

    N = div(length(ks), 2)
    UL =  [ LR[k] > 0 ? sqrt( 2/(N + 1)) * sin( pi * ks[k] * j/(N+ 1)) : 0 for k in eachindex(ks)]
    UR =  [ LR[k] < 0 ? sqrt( 2/(N + 1)) * sin( pi * ks[k] * j/(N+ 1)) : 0 for k in eachindex(ks)]

    return UL, UR
end 

function add_specific_int!(sys:: DPT_mixed, res)

    lower = get_systotal(sys) - 1
    # couple range L
    for j in 1:couple_range(sys)

        UL, UR = Uk(j, sys)
        UNN = UL .* UL' + UR .* UR'

        for k = 1:L(sys) + R(sys)
            for l =1:L(sys) + R(sys)

                res += U(sys) * UNN[k, l], "N", lower, "Cdag", k, "C", l
                #if LR(sys)[k] == LR(sys)[l]

                    #res += U(sys) * Uk( ks(sys)[k], j, sys) * Uk( ks(sys)[l], j, sys), "N", lower, "Cdag", k, "C", l
                    #res += - 1/2 * Uk( ks(sys)[k], j, sys) * Uk( ks(sys)[l], j, sys), "Cdag", k, "C", l

                    
                #end 

            end 
        end 
    end 

    return res
end 