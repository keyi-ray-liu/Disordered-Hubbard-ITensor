

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


    if includeU(sys)

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
        # if only mixed include these sites

        for j in 1:couple_range(sys)

            UL, UR = Uk(j, sys)
            UNN = UL .* UL' + UR .* UR'


            if U(sys) != 0
                for k = 1:L(sys) + R(sys)
                    for l =1:L(sys) + R(sys)

                        if UNN[k, l] != 0 
                            res += U(sys) * UNN[k, l], "N", dd_lower(sys), "Cdag", k, "C", l
                            res -= 1/2 *U(sys) * UNN[k, l],  "Cdag", k, "C", l
                        end 

                    end 
                end 
            end 
        end 

    else

        # connection to coupled site
        UL, UR = Uk(couple_range(sys), sys)

        # we separate the cases as L, R true indices can be offset due to the presence of the middle section
        for k in union(1:L(sys) - couple_range(sys) , L(sys) + couple_range(sys) + 1 : L(sys) + R(sys))
            
            k_mix = k > L(sys) ? k - LR_site_offset(sys) : k

            if UL[k_mix] != 0
                res += UL[k_mix] * t_reservoir(sys), "Cdag", k, "C", L(sys) - couple_range(sys)
                res += UL[k_mix] * t_reservoir(sys), "Cdag", L(sys) - couple_range(sys), "C", k

            elseif UR[k_mix] != 0
                res += UR[k_mix] * t_reservoir(sys), "Cdag", k, "C", L(sys) + couple_range(sys)
                res += UR[k_mix] * t_reservoir(sys), "Cdag", L(sys) + couple_range(sys), "C", k
            end 

        end 



    end 


    return res
end 