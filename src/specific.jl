

Uk(j::Int, sys::DPT_mixed) = Uk(j, ks(sys), LR(sys))

"""generate the  transformation LR vectors for the reservoir(s) SPATIAL site j, with given ordering of the momentum space, with LR cross terms being zero"""
function Uk(j::Int, ks, LR)

    N = div(length(ks), 2)
    UL =  [ LR[k] > 0 ? sqrt( 2/(N + 1)) * sin( pi * ks[k] * j/(N+ 1)) : 0 for k in eachindex(ks)]
    UR =  [ LR[k] < 0 ? sqrt( 2/(N + 1)) * sin( pi * ks[k] * j/(N+ 1)) : 0 for k in eachindex(ks)]

    return UL, UR
end 

function add_specific_int!(sys:: DPT_mixed, res)

    println("Adding mixed int")
    lower = get_systotal(sys) - 1

    
    #LR connection 
    UL, UR = Uk(1, sys)
    ULR = UL .* UR'

    for k = 1:L(sys) + R(sys)
        for l =1:L(sys) + R(sys)

            if ULR[k, l] != 0 
                res += ULR[k, l], "Cdag", k, "C", l
                res += ULR[k, l], "Cdag", l, "C", k
            end 

        end 
    end 


    # couple range L
    for j in 1:couple_range(sys)

        UL, UR = Uk(j, sys)
        UNN = UL .* UL' + UR .* UR'


        if U(sys) != 0
            for k = 1:L(sys) + R(sys)
                for l =1:L(sys) + R(sys)

                    if UNN[k, l] != 0 
                        res += U(sys) * UNN[k, l], "N", lower, "Cdag", k, "C", l
                        res -= 1/2 *U(sys) * UNN[k, l],  "Cdag", k, "C", l
                    end 

                end 
            end 
        end 
    end 

    return res
end 