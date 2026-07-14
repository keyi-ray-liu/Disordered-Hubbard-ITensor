function connected_sites(nsq::Int, i::Int)
    # 1) collect all corners
    coordset = Set{Tuple{Int,Int}}()
    for k in 0:nsq-1
        base = k % 2
        for dx in (0,1), dy in (0,1)
            push!(coordset, (k + dx, base + dy))
        end
    end

    # 2) sort row‐major: first by y, then by x
    coords = sort(collect(coordset), by = c -> (c[2], c[1]))

    # 3) bounds‐check
    N = length(coords)
    @assert 1 ≤ i ≤ N "site index i=$i out of bounds 1:$N"

    # 4) brute‐force neighbor scan
    x,y = coords[i]
    nbrs = Int[]
    for (j,(x2,y2)) in enumerate(coords)
        if (abs(x2-x)==1 && y2==y) || (abs(y2-y)==1 && x2==x)
            push!(nbrs, j)
        end
    end

    return sort(nbrs)
end


function connected_sites2(nsq::Int, i::Int)
    # 1) collect all corner coords in a Set to dedupe
    coordset = Set{Tuple{Int,Int}}()
    for k in 0:nsq-1
        base = k % 2
        for dx in (0,1), dy in (0,1)
            push!(coordset, (k + dx, base + dy))
        end
    end

    # 2) turn into a sorted Vector of coords (row-major: y first, then x)
    coords = sort(collect(coordset), by = c -> (c[2], c[1]))

    # 3) build the list of edges from each square’s four sides
    edges = Tuple{Int,Int}[]
    for k in 0:nsq-1
        base = k % 2
        # the four corners of square k
        a = (k,   base)
        b = (k+1, base)
        c = (k,   base+1)
        d = (k+1, base+1)
        # its bonds: (a–b), (a–c), (b–d), (c–d)
        for (p,q) in ((a,b),(a,c),(b,d),(c,d))
            # find their 1-based indices by linear scan
            ip = findfirst(x->x==p, coords)
            iq = findfirst(x->x==q, coords)
            # store as an (min,max) tuple to avoid directionality
            push!(edges, (min(ip,iq), max(ip,iq)))
        end
    end
    # remove duplicate edges
    edges = unique(edges)

    # 4) sanity‐check site index
    N = length(coords)
    @assert 1 ≤ i ≤ N "site index i=$i out of range 1:$N"

    # 5) collect neighbors by looking at all edges
    nbrs = Int[]
    for (u,v) in edges
        if u == i
            push!(nbrs, v)
        elseif v == i
            push!(nbrs, u)
        end
    end

    return sort(nbrs)
end

coords = connected_sites2(4, 2)
println(coords)