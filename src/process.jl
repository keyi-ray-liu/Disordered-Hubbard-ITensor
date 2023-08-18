"""Get number of maximum sweeps"""
function getmaxsweep(;L=12, tol=1e-11)
  maxsweepdim = 500

  output = @capture_out begin
    paras = setpara(L=L, N=Int(L//2), sweepdim=maxsweepdim)
    main(paras;)
  end

  maxsweep = stringmaxproc(output, tol)
  return maxsweep
end

"""work function that get max sweeps for different lengths of chains"""
function sweepstat()
  res = zeros(0)
  for L in range(12, 100, step=2)
    append!(res, getmaxsweep(L=L))
    println(res)
  end
end

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

function cal_overlap()

  workdir = getworkdir()
  baseline = h5open(workdir * "QE.h5", "r")
  energy = readdlm( workdir * "energy")


  start = 0.1
  fin = 20.0
  step = 0.1

  states = sort(keys(baseline), by= x -> parse(Int, x[4:end]))

  print(states)
  ex = length(states)

  if ex != length(energy)
    error("# of energies do not match # of eigenstates")
  end 

  exwf = Vector{MPS}(undef, ex)

  for (i, key) in enumerate(states)
    exwf[i] = read(baseline, key, MPS)
  end 

  close(baseline)

  L = start:step:fin

  overlap = Array{ComplexF64}(undef, (length(L), ex))
  overlapnorm = zeros((length(L), ex))

  for (i, t) in enumerate(L)

    TE = h5open(workdir * "t" * string(t) *".h5", "r")
    ψ = read( TE, "psi", MPS)
    close(TE)

    for (j, wf) in enumerate(exwf)

      cur = inner(ψ', wf)
      cur *= exp( -im * energy[j] * (t + step))

      overlap[i, j] = cur
      overlapnorm[i, j] = abs(cur)
      println(i, j, " done")

    end 
  end 

  writedlm( workdir * "overlap", overlap)
  writedlm( workdir * "overlapnorm", overlapnorm)
  
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
function temp_occ(num)
  
  
  function get_time(raw::String)

    rawt = SubString(raw, 1 + 1 + length(method), length(raw) - 3)
    t = parse(Float64, rawt )

    return t
  end 

  workdir = getworkdir()
  control = parse(Int, num)


  if control == 1
    method = "TEBD"

  elseif control == 2
    method = "TDVP"

  end 

  res = []
  bonds = []
  T = []


  files = filter(x->occursin(method,x), readdir(workdir))
  sort!(files, by=get_time)

  for file in files

    wf = h5open(workdir * file, "r")
    ψ = read( wf, "psi", MPS)
    close(wf)

    t = get_time(file)
    bond = checkmaxbond(ψ)
    occ = expect(ψ, "N")
      
    append!(T, t)
    append!(bonds, bond)
    append!(res, [occ])
  end 

  writedlm(workdir *"time", T)
  writedlm(workdir * "occ", res)
  writedlm(workdir * "bonddim", bonds)

  gs = h5open( workdir * "gs.h5", "r")
  gs_wf = read(gs, "psi1", MPS)
  close(gs)

  occ_gs = expect(gs_wf, "N")
  writedlm( workdir * "gs", occ_gs)
end 

function eigen_overlap()

  workdir = getworkdir()
  gs = h5open(workdir * "target.h5", "r")
  gswf = read( gs, "psi1", MPS)

  print(length(gswf))
  total = 0 
  wf = h5open(workdir * "wf.h5", "r")
  exst = sort(keys(wf), by= x-> parse(Int, x[4:end]))
  for ex in eachindex(exst)
    psi = read(wf, exst[ex], MPS )


    overlap = inner(psi', gswf)
    total += abs(overlap) ^2
    println(total)
  end 

end 


function get_eigen_occ(tcd_dict, phase_overlap)

  num = minimum( map(maximum, collect(zip(collect(keys(tcd_dict))...))))

  println(num)
  res = zeros( length(phase_overlap[1]))

  for left in 1:num
    for right in left+1:num

      tcd = tcd_dict[ (left, right)]
      res = res .+ conj(phase_overlap[left]) .* tcd .* phase_overlap[right]
      res = res .+ conj(phase_overlap[right]) .* conj.(tcd) .* phase_overlap[left]

    end 
  end 
  
  for idx in 1:num
    tcd = tcd_dict[ (idx, idx)]
    res = res.+ abs2( phase_overlap[idx]) * tcd
  end 

  return real.(res)

end 
