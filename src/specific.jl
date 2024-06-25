
Uk(j::Int, sys::DPT_avg) = Uk(j, ks(sys), LR(sys))
Uk(j::Int, sys::DPT_mixed) = Uk(j, ks(sys), LR(sys))
Uk(j::Int, sys::DPT_graph) = typeof(sys.dpt) == DPT_mixed ? Uk(j, sys.dpt) : error("subsystem not mixed")


"""generate the  transformation LR vectors for the reservoir(s) SPATIAL site j, with given ordering of the momentum space, with LR cross terms being zero"""
function Uk(j::Int, ks, LR)

    N = div(length(ks), 2)
    UL =  [ LR[k] > 0 ? sqrt( 2/(N + 1)) * sin( pi * ks[k] * j/(N+ 1)) : 0 for k in eachindex(ks)]
    UR =  [ LR[k] < 0 ? sqrt( 2/(N + 1)) * sin( pi * ks[k] * j/(N+ 1)) : 0 for k in eachindex(ks)]

    return UL, UR
end 

function add_specific_int!(sys:: Union{DPT_mixed, DPT_graph}, res)


    if typeof(sys) == DPT_graph && typeof(sys.dpt) == DPT
        return res
    end 

    @info "Adding mixed int"

    if includeU(sys)

        #LR connection 
        UL, UR = Uk(1, sys)
        ULR = UL .* UR'

        for (k_mix, k_actual) in enumerate(union(L_begin(sys): L_end(sys) , R_begin(sys) : R_end(sys)))
            for (l_mix, l_actual) in enumerate(union(L_begin(sys): L_end(sys) , R_begin(sys) : R_end(sys)))

                if ULR[k_mix, l_mix] != 0 
                    res += ULR[k_mix, l_mix], "Cdag", sitemap(sys, k_actual), "C", sitemap(sys, l_actual)
                    res += ULR[k_mix, l_mix], "Cdag", sitemap(sys, l_actual), "C", sitemap(sys, k_actual)
                end 

            end 
        end 


        # couple range L
        # if only mixed include these sites

        for j in 1:couple_range(sys)

            UL, UR = Uk(j, sys)
            UNN = UL .* UL' + UR .* UR'


            if U(sys) != 0
                for (k_mix, k_actual) in enumerate(union(L_begin(sys): L_end(sys) , R_begin(sys) : R_end(sys)))
                    for (l_mix, l_actual) in enumerate(union(L_begin(sys): L_end(sys) , R_begin(sys) : R_end(sys)))

                        if UNN[k_mix, l_mix] != 0 
                            res += U(sys) * UNN[k_mix, l_mix], "N", sitemap(sys, dd_lower(sys)), "Cdag", sitemap(sys, k_actual), "C", sitemap(sys, l_actual)
                            res -= 1/2 *U(sys) * UNN[k_mix, l_mix],  "Cdag", sitemap(sys, k_actual), "C", sitemap(sys, l_actual)
                        end 

                    end 
                end 
            end 
        end 

    else

        # connection to coupled site
        #UL, UR = Uk(couple_range(sys), sys)
        UL, UR = Uk(1, sys)

        # we separate the cases as L, R true indices can be offset due to the presence of the middle section
        for (k_mix, k_actual) in enumerate(union(L_begin(sys): L_contact(sys) - 1 , R_contact(sys) + 1 : R_end(sys)))
            
            if UL[k_mix] != 0
                res += UL[k_mix] * t_reservoir(sys), "Cdag", sitemap(sys, k_actual), "C", sitemap(sys, L_contact(sys))
                res += UL[k_mix] * t_reservoir(sys), "Cdag", sitemap(sys, L_contact(sys)), "C", sitemap(sys, k_actual)

            elseif UR[k_mix] != 0
                res += UR[k_mix] * t_reservoir(sys), "Cdag", sitemap(sys, k_actual), "C", sitemap(sys, R_contact(sys))
                res += UR[k_mix] * t_reservoir(sys), "Cdag", sitemap(sys, R_contact(sys)), "C", sitemap(sys, k_actual)
            end 

        end 



    end 


    return res
end 



function add_specific_int!(sys:: DPT_avg, res)

    # spatial
    if typeof(sys.dpt) == DPT_mixed

        @info "Adding DPTavg mixed int"


        # connection to coupled site
        #UL, UR = Uk(couple_range(sys), sys)
        UL, UR = Uk(1, sys)

        # we separate the cases as L, R true indices can be offset due to the presence of the middle section
        for (k_mix, k_actual) in enumerate(union(L_begin(sys): L_contact(sys) - 1 , R_contact(sys) + 1 : R_end(sys)))
            
            if UL[k_mix] != 0
                res += UL[k_mix] * t_reservoir(sys), "Cdagup", sitemap(sys, k_actual), "Cup", sitemap(sys, L_contact(sys))
                res += UL[k_mix] * t_reservoir(sys), "Cdagup", sitemap(sys, L_contact(sys)), "Cup", sitemap(sys, k_actual)

            elseif UR[k_mix] != 0
                res += UR[k_mix] * t_reservoir(sys), "Cdagup", sitemap(sys, k_actual), "Cup", sitemap(sys, R_contact(sys))
                res += UR[k_mix] * t_reservoir(sys), "Cdagup", sitemap(sys, R_contact(sys)), "Cup", sitemap(sys, k_actual)
            end 

        end 

    end 


    return res
end 