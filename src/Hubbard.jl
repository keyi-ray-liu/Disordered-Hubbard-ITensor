HubbardRepulsion(sys::Systems, j::Int) = systype(sys) == "Fermion" ? 0.0 : U(sys)

function HubbardRepulsion(sys::SD_array, j::Int)

    source = get_systotal(sys.source)
    array = get_systotal(sys.array)

    if source < j <= source + array
        return U(sys.array)

    else
        return 0
    end 

end 

function add_HubbardRepulsion!(sys::Systems, res::OpSum)
    
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
