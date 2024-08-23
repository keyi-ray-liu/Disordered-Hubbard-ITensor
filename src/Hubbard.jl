HubbardRepulsion(sys::systems, j::Int) = 0.0
HubbardRepulsion(sys::NF_square, j::Int) = U(sys)


function add_HubbardRepulsion!(sys::systems, res::OpSum)
    
    @info "Adding HubbardRepulsion"
    systotal = get_systotal(sys)

    sys_type = systype(sys)

    if sys_type == "Fermion"
        op = "I"

    elseif sys_type == "Electron"
        op = "Nupdn"

    end 

    for j=1 :systotal
        # E-E and N-
        if HubbardRepulsion(sys, j) != 0
            res += HubbardRepulsion(sys, j), op, sitemap(sys, j)
        end 
    end 

    return res

end 
