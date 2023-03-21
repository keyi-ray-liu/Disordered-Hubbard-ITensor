"""worker function perfomring time evolution"""
function time_evolve(ψ, sites, paras, start, fin, τ)

    QE = paras["QE"]
    ζ_ne, ζ_ee = paras["ζ"]
    L = paras["L"]
    λ_ne = paras["int_ne"]
    λ_ee = paras["int_ee"]
    CN = paras["CN"]
    t = paras["t"]
    decay = paras["decay"]
    QEoffset = paras["QEoffset"]
    QN = paras["QN"]
    exch = paras["exch"]
    dp = paras["dp"]
    ζ_dp = paras["ζ_dp"]
    disorder = paras["disorder"]
    range = paras["range"]
    self = paras["self_nuc"]
    QEen = paras["QEen"]
    cutoff = paras["TEcutoff"]
    maxdim = paras["TEdim"]
    method = para["TEmethod"]

    if QE == 0
        throw(ArgumentError("Dynamics have to include quantum emitters"))

    elseif QE > 2
        throw(ArgumentError("More than two QE, not allowed"))
    end 

    head = (QE > 0) * (QN + 1)
    # get disorder (if no then zeros) 
    disx, disy = setdisorder(disorder, L)

    gates = ITensor[]
    
    λ_ne = λ_ne * CN / L

    # two-site hopping gates

    # TEBD approximation factor
    factor = 2.0
    for j in 1:(L - 1)

        r = dis(j, j + 1, disx, disy)
        hop = hopping(decay, r)

        s1 = sites[j + head ]
        s2 = sites[j + head + 1]
        hj =
          - t * hop * op("C", s1) * op("Cdag", s2) +
          - t * hop * op("C", s2) * op("Cdag", s1)

        gatefy!(gates, factor, hj, τ)
    end

    # two site NN gates
    for j=1 :L 

        s1 = sites[j + head]
        # E-E
        for k= max(1, j - range): j-1
            
            s2 = sites[k + head]
            # delta function setting up the exchange between nearest neighbor
            ifexch = (1 -  ==(1, abs(j - k)) * exch )
            
            # add the ee interaction term one by one
            eejk = λ_ee * ifexch / ( dis(j, k, disx, disy) + ζ_ee) * op("N", s1) * op("N", s2)

            gatefy!(gates, factor, eejk, τ)

        end

    end 

    # add three site dp gates:
    if QE > 0

        # site 1. in QN: this is the AUX site. otherwise this is the QE site
        s0 = sites[1]
        dp_left = dp[1]
        ζ_dp_left = ζ_dp[1]
        cavg = CN / L 

        for left = 1 : L 
            
            s2 = sites[left + head]

            r = dis(left, QEoffset, disx, disy)
            r_dp = r ^ 3 + ζ_dp_left

            if !QN
                dpla = dp_left * r  / r_dp * op( "x", s0) * op("N", s2)
                dpl = -dp_left * r *  cavg / r_dp * op( "x", s0)

            else
                s1 = sites[head]

                #two level with AUX site
                dpla = dp_left * r / r_dp  * op( "C", s0) * op("Cdag", s1) * op("N", s2) + 
                dp_left * r / r_dp  * op( "C", s1) * op("Cdag", s0) * op("N", s2)
    
                # the offset term, to set the 'center of mass'
                # calculated as a uniform distribution:  L * N / 2
                dpl =  - dp_left * r * cavg / r_dp *  op("C", s0) * op( "Cdag", s1) - 
                dp_left * r * cavg / r_dp *  op("C", s1) * op( "Cdag", s0) 

            end 

            gatefy!(gates, factor, dpla, τ)
            gatefy!(gates, factor, dpl, τ)
          end 

    end 

    if QE > 1
        # in this case, s0 is the rightmost site, either L + 2 or L + 4
        s0 = sites[L + head * 2]
        dp_right = dp[1]
        ζ_dp_right = ζ_dp[1]
        
        for right = 1 : L 
            
            s2 = sites[right + head]

            r = ( L - 1 + QE * (QEoffset + 1) -  dis(right, QEoffset, disx, disy))  
            r_dp = r ^ 3 + ζ_dp_right

            if !QN
                dpra = dp_right * r  / r_dp * op( "x", s0) * op("N", s2)
                dpr = -dp_right * r *  cavg / r_dp * op( "x", s0)

            else
                s1 = sites[L + head + 1]

                #two level with AUX site
                dpra = dp_right * r / r_dp  * op( "C", s0) * op("Cdag", s1) * op("N", s2) + 
                dp_right * r / r_dp  * op( "C", s1) * op("Cdag", s0) * op("N", s2)
    
                # the offset term, to set the 'center of mass'
                # calculated as a uniform distribution:  L * N / 2
                dpr =  - dp_right * r * cavg / r_dp *  op("C", s0) * op( "Cdag", s1) - 
                dp_right * r * cavg / r_dp *  op("C", s1) * op( "Cdag", s0) 

            end 

            gatefy!(gates, factor, dpra, τ)
            gatefy!(gates, factor, dpr, τ)
          end 

    end 


    # reverse gates
    append!(gates, reverse(gates))

    # two(three)-site gate end
    # add one-site TE gates

    factor = 1.0
    for j=1 :L 

        s1 = sites[j + head]
        
        # N-E
    
        # check if self self_nuc
    
        if self
            cursum = 0
        else
            cursum = -λ_ne / ζ_ne
        end
    
        for l=1 :L 
            cursum += λ_ne / ( dis(j, l, disx, disy) + ζ_ne)
        end
    
        nej = cursum * op("N", s1)
        gatefy!(gates, factor, nej, τ)
    end

    # QE diagonal energy
    if QE > 0
        g = op("N", sites[1]) * QEen
        gatefy!(gates, factor, g, τ)
    end 
    
    if QE > 1
        g = op("N", sites[L + 2]) * QEen
        gatefy!(gates, factor, g, τ)
    end 
    
    # add QN protecting gates 
    if !QN
        Λ = 30.0
    
        for i= 1:L

            s1 = sites[i + head]
        
            # linear gate
            li = - 2 * Λ * N * op( "N", s1)
            gatefy!(gates, factor, li, τ)

            # quadratic gate
            for j =1:L
                
                s2 = sites[j + head]

                # temp fix: issue with double gates
                if s1 != s2
                    lij = Λ * op("N", s1) * op("N", s2)
                end 

                gatefy!(gates, factor, lij, τ)

            end 
    
        end 
    
    end 

    # one site gate end

    prefix = getworkdir()


    for dt in start:τ:fin
        #Sz = expect(psi, "Sz"; sites=c)

        # testing the occupation on the whole chain plus emitter
        println("$dt")
    
        ψ = apply(gates, ψ; cutoff=cutoff, maxdim = maxdim)
        normalize!(ψ)

        wf = h5open( prefix * "t" * string(dt) * ".h5", "w")
        write(wf, "psi", ψ)
        close(wf)
        
    end


    return nothing

end 