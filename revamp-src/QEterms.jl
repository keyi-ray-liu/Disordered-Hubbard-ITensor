
"""default no QE"""
QECoupling(sys::systems, j) = []


function QECoupling(sys::QE_two, j) 

    # no coupling for QE sites themselves, as they are not coupled to other QEs
    if j <= QESITES || j > get_systotal(sys) - QESITES
        return []

    else

        r_qe1 = dis(j, QESITES + 1, sys)
        qe1 = dp(sys) * r_qe1/ ( r_qe1^3 + ζ(sys) )

        r_qe2 = dis(j, get_systotal(sys) - QESITES, sys)
        qe2 = dp(sys) * r_qe2/ ( r_qe2^3 + ζ(sys))

        return [ [qe1, 1], [qe2, get_systotal(sys) - 1]]

    end 
    

end 

"""For the flat X QE, currently chain sites only couple to their respective QE"""
function QECoupling(sys::QE_flat_SIAM, j) 

    # no coupling for QE sites themselves, as they are not coupled to other QEs
    qe_loc, chain_begin, chain_end = get_sys_loc(sys, j)

    if chain_begin <= j <= chain_end

        if j > qe_loc 
            r_qe = dis(j, chain_begin, sys)
        else
            r_qe = dis(j, chain_end, sys)
        end 

        qe = dp(sys)  * r_qe / (r_qe^3 + ζ(sys))
        res = [[qe, qe_loc]]

    else
        res = []

    end 

    return res

end 


"""
All logic are wrapped in respective functions. The QE diagonal Energies are wrapped in onsite function, where the offset QE hopping should be included in the hopping part of the Hamiltoian
"""
function add_qe!(sys::systems, res::OpSum)

    println("Adding QE")

    sys_type = type(sys)
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

                res += dp, op1, kex, op2, kgs,  opn, j
                res += dp, op1, kgs, op2, kex, opn, j

                res += - dp * offset_scale(sys), op1, kex, op2, kgs
                res += - dp * offset_scale(sys), op1, kgs, op2, kex
            end 
        end 

    end 

    return res

end



