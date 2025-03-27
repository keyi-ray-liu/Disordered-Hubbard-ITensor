
"""default no QE"""
QECoupling(sys::Systems, j) = []


function QECoupling(sys::QE_two, j) 

    # no coupling for QE sites themselves, as they are not coupled to other QEs
    if j <= QESITES || j > get_systotal(sys) - QESITES
        return []

    else

        r_qe1 = qedis(j, QESITES + 1, sys)
        qe1 = dp(sys) * r_qe1/ ( r_qe1^3 + ζ(sys) )

        r_qe2 = qedis(j, get_systotal(sys) - QESITES, sys)
        qe2 = dp(sys) * r_qe2/ ( r_qe2^3 + ζ(sys))

        return [ [qe1, 1], [qe2, get_systotal(sys) - 1]]

    end 
    

end 

# function QECoupling(sys::QE_parallel, j)

#     uppertotal = get_uppertotal(sys)
#     systotal = get_systotal(sys)
#     center_lower = systotal - 1

#     if j <= uppertotal
#         qe =  QECoupling(sys.upper, j)

#     elseif j < center_lower
#         qe =  [ [U, k + uppertotal] for (U, k) in QECoupling(sys.lower, j - uppertotal)]

#     # no QE on contact
#     else
#         qe = []
#     end 

#     return qe

# end 


function QECoupling(sys::QE_HOM, j)

    uppertotal = get_uppertotal(sys)
    systotal = get_systotal(sys)

    if j <= uppertotal
        qe =  QECoupling(sys.upper, j)

    elseif j < systotal
        qe =  [ [U, k + uppertotal] for (U, k) in QECoupling(sys.lower, j - uppertotal)]

    else
        qe = [] 
    end 

    return qe

end 


# """For the flat X QE, currently chain sites only couple to their respective QE"""
# function QECoupling(sys::QE_flat_SIAM, j) 

#     # no coupling for QE sites themselves, as they are not coupled to other QEs
#     qe_loc, chain_begin, chain_end = get_sys_loc(sys, j)

#     if chain_begin <= j <= chain_end

#         if j > qe_loc 
#             r_qe = qedis(j, chain_begin, sys)
#         else
#             r_qe = qedis(j, chain_end, sys)
#         end 

#         qe = dp(sys)  * r_qe / (r_qe^3 + ζ(sys))
#         res = [[qe, qe_loc]]

#     else
#         res = []

#     end 

#     return res

# end 

# QECoupling(sys::QE_G_SIAM, j) = QECoupling(sys.system, j)


"""
All logic are wrapped in respective functions. The QE diagonal Energies are wrapped in onsite function, where the left_offset QE hopping should be included in the hopping part of the Hamiltoian
"""
function add_qe!(sys::Systems, res::OpSum)

    @info "Adding QE"

    sys_type = systype(sys)
    systotal = get_systotal(sys)
    

    if sys_type == "Electron"
        opn = "Ntot"
        opc = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]
    
    elseif sys_type == "Fermion"
        opn = "N"
        opc = [ ["C", "Cdag"]]

    end 

    for j in 1:systotal

        for (dp, kex) in QECoupling(sys, j)
            
            kex = trunc(Int, kex)
            kgs = kex + 1
            for op in opc
                op1, op2 = op

                res += dp, op1, sitemap(sys, kex), op2, sitemap(sys, kgs),  opn, sitemap(sys, j)
                res += dp, op1, sitemap(sys, kgs), op2, sitemap(sys, kex), opn, sitemap(sys, j)

                res += - dp * offset_scale(sys), op1, sitemap(sys, kex), op2, sitemap(sys, kgs)
                res += - dp * offset_scale(sys), op1, sitemap(sys, kgs), op2, sitemap(sys, kex)
            end 
        end 

    end 

    return res

end



