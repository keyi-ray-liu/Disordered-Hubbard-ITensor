function HoppingNeighbor(sys::QE_flat_SIAM, j)
end 


"""
Here, we have either the offset hopping in the QE, or the chain hopping
"""
function HoppingNeighbor(sys::QE_two, j)

    systotal = get_systotal(sys)
    t = t(sys)

    # we only hop onwards, no hop in QE and last site
    if j <3 || j > systotal - 3
        hop = []
    else
        hop = [[t, j + 1]]
    end 

    return hop

end 

# for a flat chain, return next NN until end
function HoppingNeighbor(sys::Chain_only, j)

    dims = dims(sys)
    t = t(sys)
    return j < dims ? [[t, j + 1] ]: []

end 