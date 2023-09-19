"""work function that performs data io and calculates exp vals of observables, in this case, TCD"""
function cal_observe(;key="wf.h5")

  workdir = getworkdir()
  obs = "TCD"

  # set range of ex-ex calculations
  cutoff1 = 50
  cutoff2 = 50

  for file in filter(x->occursin(key,x), readdir(workdir))

    wf = h5open( workdir * file, "r")
    exst = sort!(keys(wf), by= x-> parse(Int, x[4:end]))
    numst = length(exst)
    psis = []

    for ex in exst
      psi = read(wf, ex, MPS)
      append!( psis, [psi])

    end

    for i = 1:min(cutoff1, numst)
      for j = i:min(cutoff2, numst)

        suffix = file[3:4] * "_ex_" * string(i) * "_" *  string(j)
        #cccorr = correlation_matrix(psi,"Cdag","C")

        if obs == "CC"
          obsval = cal_corr(psis[i], psis[j])

        elseif obs == "TCD"
          obsval = cal_tcd(psis[i], psis[j])

        else
          println("unknown observable, abort")
          exit()
        end 

        writedlm( workdir * obs * "_" *  suffix, obsval)
      end
    end
  end

end


function gs_occ()

    workdir = getworkdir()
  
    gs = h5open( workdir * "target.h5", "r")
    gs_wf = read(gs, "psi1", MPS)
    close(gs)
  
    occ_gs = expect(gs_wf, "N")
    writedlm( workdir * "gs", occ_gs)
end 
  

"""calculates the occ for the available t slices"""
function time_obs(para)
  
  
  function get_time(raw::String)

    rawt = SubString(raw, 1 + 1 + length(method), length(raw) - 3)
    t = parse(Float64, rawt )

    return t
  end 

  workdir = getworkdir()
  method = para["method"]
  obs = para["obs"]

  if obs == "occ"
    occs = []
    bonds = []

  elseif obs == "tcd"

    qewf, _, _ = load_qe()
    λ_ee = para["λ_ee"]
    ζ = para["ζ"]

  elseif obs == "current"
    
    current = []

  else
    error("Unrecognized obs")
  end 

  T = []


  files = filter(x->occursin(method,x), readdir(workdir))
  sort!(files, by=get_time)

  for file in files

    wf = h5open(workdir * file, "r")
    ψ = read( wf, "psi", MPS)
    close(wf)

    t = get_time(file)
    append!(T, t)

    if obs == "occ"
      bond = checkmaxbond(ψ)
      occ = expect(ψ, "N")
      append!(bonds, bond)
      append!(occs, [occ])

    elseif obs == "tcd"

        if !isfile( workdir * "TCD_te_RE$t")

            tcds = []
            println("cal tcd $t")

            for ex in eachindex(qewf)
            append!( tcds, [cal_tcd( qewf[ex], ψ)])
            end 
        
        else

            tcd_re = readdlm( workdir * "TCD_te_RE$t" )
            tcd_im = readdlm( workdir * "TCD_te_IM$t")

            tcds = tcd_re + 1im * tcd_im

        end 

        gpi = get_gpi(tcds, λ_ee, ζ)

        writedlm( workdir * "gpi_RE$t", real(gpi))
        writedlm( workdir * "gpi_IM$t", imag(gpi))
        writedlm( workdir * "TCD_te_RE$t", real(tcds))
        writedlm( workdir * "TCD_te_IM$t", imag(tcds))


    elseif obs == "current"

        I = cal_current(ψ, para)
        println("current at $t", I)
        append!(current, I)

    end 
    
  end 

  writedlm(workdir *"time", T)

  if obs == "occ"
    writedlm(workdir * "occ", occs)
    writedlm(workdir * "bonddim", bonds)

    gs_occ()

  elseif obs == "current"

    writedlm(workdir * "current", current)
  end


end 
