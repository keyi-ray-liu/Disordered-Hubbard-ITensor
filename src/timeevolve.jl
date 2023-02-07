"""worker function perfomring time evolution"""
function time_evolve(ψ; t=1.0, L=12, N=6, CN=6, λ_ne=1.0, λ_ee=1.0, QE=2, QN=false, QEen=0.4, dp=0.25, time=5, τ=0.1, disorder=false, range=1000, self=false, cutoff=1E-8)

    if QE == 0
        throw(ArgumentError("Dynamics have to include quantum emitters"))

    elseif QE > 2
        throw(ArgumentError("More than two QE, not allowed"))
    end 

    head = (QE > 0)

    # get disorder (if no then zeros) 
    disx, disy = setdisorder(disorder, L)

    sites =  siteinds("Fermion", N + QE; conserve_qns=true)
    gates = ITensor[]
    

    # two-site hopping gates

    factor = 2 
    for j in 1:(L - 1)

        r = dis(j, j + 1, disx, disy)
        hop = hopping(decay, r)

        s1 = sites[j + head ]
        s2 = siste[j + head + 1]
        hj =
          - t * hop * op("C", s1) * op("Cdag", s2) +
          - t * hop * op("C", s2) * op("Cdag", s1)

        gatefy(gates, factor, hj, τ)
    end

    # add three site dp gates:
    if QE > 0

        for left = 1 : L 
            
            s0 = sites[1]
            s1 = sites[left + head]

            r = dis(left, disx, disy)
            for all = 1 : L
                
                s2 = sites[all + head]
                # r0 determines the overall 'weight' of sites
                r0 = dis(all, disx, disy)
                
                #off-diagonal two level transition term
                dpla = - dp * r0 / r ^ 3 * op( "x", s0) * op("N", s1) * op("N", s2)
                gatefy(gates, factor, dpla, τ)
              
            end 
            
            # the offset term, to set the 'center of mass'
            # calculated as a uniform distribution:  L * N / 2
            dpl =  dp * L * N / ( 2 * r^3) *  op("x", s0) * op( "N", s1)
            gatefy(gates, factor, dpl, τ)
          end 

    end 

    if QE > 1

        for right = 1 : L 
            
            s0 = sites[L + 2]
            s1 = sites[right + head]

            r = dis(right, disx, disy)
            for all = 1 : L
                
                s2 = sites[all + head]
                # r0 determines the overall 'weight' of sites
                r0 = dis(all, disx, disy)
                
                #off-diagonal two level transition term
                dpla = - dp * (L + 1 - r0) / ( L + 1 - r) ^ 3 * op( "x", s0) * op("N", s1) * op("N", s2)
                gatefy(gates, factor, dpla, τ)
              
            end 
            
            # the offset term, to set the 'center of mass'
            # calculated as a uniform distribution:  L * N / 2
            dpl =  dp * L * N / ( 2 * r^3) *  op("x", s0) * op( "N", s1)
            gatefy(gates, factor, dpl, τ)
          end 

    end 


    # reverse gates
    append!(gates, reverse(gates))
    # two-site gate end



    # add one-site TE gates

    factor = 1
    for j=1 :L 

        s1 = sites[j + head]
        # E-E
        for k= max(1, j - range):j-1
            
            s2 = sites[k + head]
            # delta function setting up the exchange between nearest neighbor
            ifexch = (1 -  ==(1, abs(j - k)) * exch )
            
            # add the ee interaction term one by one
            eejk = λ_ee * ifexch / ( dis(j, k, disx, disy) + ζ_ee) * op("N", s1) * op("N", s2)

            gatefy(gates, factor, eejk, τ)

        end
        
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
        gatefy(gates, factor, nej, τ)
    end

    # QE diagonal energy
    if QE > 0
        g = op("N", sites[1]) * QEen
        gatefy(gates, factor, g, τ)
    end 
    
    if QE > 1
        g = op("N", sites[L + 2]) * QEen
        gatefy(gates, factor, g, τ)
    end 
    
    # add QN protecting gates 
    if !QN
        Λ = 30.0
    
        for i= 1:L

            s1 = sites[i + head]
        
            # linear gate
            li = - 2 * Λ * N * op( "N", s1)
            gatefy(gates, factor, li, τ)

            # quadratic gate
            for j =1:L
                
                s2 = sites[j + head]
                lii = Λ * op("N", s1) * op("N", s2)
                gatefy(gates, factor, lii, τ)

            end 
    
        end 
    
    end 

    # one site gate end

    # further preparation of the initial state, if needed
    ψ = TE_stateprep(ψ, QE)

    for dt in 0.0:τ:time
        #Sz = expect(psi, "Sz"; sites=c)
        println("$dt")
        dt = time && break
    
        ψ = apply(gates, ψ; cutoff)
        normalize!(ψ)
    end



end 