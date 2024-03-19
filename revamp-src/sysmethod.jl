

function gen_state(sys::systems)

    state_str = gen_state_str(sys)


    

end 

function gen_hamiltonian(sys::systems)

    res = OpSum()

    res = add_hop(sys, res)
    res = add_DensityDensity(sys, res)
    res = add_onsite(sys, res)
    res = add_qe(sys, res)
    res = add_HubbardRepulsion(sys, res)
    res = add_specific_int(sys, res)

    return res
end 


# All the system specific logic is wrapped with multiple-dispatch neighbor functions, we no longer perform these logics within the code
function add_DensityDensity(sys::systems, res::OpSum)
    

    println("Adding EE and NE")
    type = type(sys)

    systotal = get_systotal(sys)

    if type == "Fermion"
        ops = "N"
  
    elseif type == "Electron"
        ops = "Ntot"
  
    end 

    for j=1 :systotal
        # E-E and N-E

        for (U, k) in DenDenNeighbor(sys, j)
            
            # delta function setting up the exchange between nearest neighbor

            # because k < j, we check if j is in the nn of k
            res += U, ops, j, ops, k

        end

    end 

    return res

end 

"""This function only concerns with the 'simple' hopping, complex terms such as QE dipole offsets are calculated elsewhere"""
function add_hop(sys::systems, res::OpSum)
    

    println("Adding all hopping")
    type = type(sys)
    systotal = get_systotal(sys)

    if type == "Fermion"
        operators = [ ["C", "Cdag"]]

    elseif type == "Electron"
        operators = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]

    end 

    for j in systotal

        for (v, k) in HoppingNeighbor(sys, j)

            for operator in operators

                op1, op2 = operator

                res += v, op1, j, op2, k
                res += v, op1, k, op2, j
            end 

        end 
    end 
    

    return res

end 


"""
All logic are wrapped in respective functions. The QE diagonal Energies are wrapped in onsite function, where the offset QE hopping should be included in the hopping part of the Hamiltoian
"""
function add_qe(sys::systems, res::OpSum)

    println("Adding QE")

    type = type(sys)
    systotal = get_systotal(sys)
    

    if type == "Electron"
        opn = "Ntot"
        opc = [ ["Cup", "Cdagup"], ["Cdn", "Cdagdn"]]
    
    elseif type == "Fermion"
        opn = "N"
        opc = [ ["C", "Cdag"]]

    end 

    for j in 1:systotal

        for (dp, kex) in QECoupling(sys, j)
            
            kgs = kex + 1
            for op in opc
                op1, op2 = op

                res += dp, op1, kex, op2, kgs,  opn, j
                res += dp, op1, kgs, op2, kex, opn, j
            end 
        end 

    end 


end


function add_onsite(sys::systems, res::OpSum)
    
    println("Adding onsite")
    type = type(sys)

    systotal = get_systotal(sys)

    if type == "Fermion"
        ops = "N"
  
    elseif type == "Electron"
        ops = "Ntot"
  
    end 

    for j=1 :systotal
        # E-E and N-
        res += Onsite(sys, j), ops, j

    end 
    return res

end 


function add_HubbardRepulsion(sys::systems, res::OpSum)
    
    println("Adding HubbardRepulsion")
    systotal = get_systotal(sys)

    for j=1 :systotal
        # E-E and N-
        res += HubbardRepulsion(sys, j), "Nupdn", j
    end 
    return res

end 
