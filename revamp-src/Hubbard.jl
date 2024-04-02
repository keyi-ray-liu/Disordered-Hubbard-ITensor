HubbardRepulsion(sys::systems, j::Int) = 0.0
HubbardRepulsion(sys::NF_square, j::Int) = U(sys)


function add_HubbardRepulsion!(sys::systems, res::OpSum)
    
    println("Adding HubbardRepulsion")
    systotal = get_systotal(sys)

    sys_type = type(sys)

    if sys_type == "Fermion"
        op = "I"

    elseif sys_type == "Electron"
        op = "Nupdn"

    end 

    for j=1 :systotal
        # E-E and N-
        res += HubbardRepulsion(sys, j), op, j
    end 

    return res

end 
