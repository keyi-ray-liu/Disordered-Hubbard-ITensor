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
  occ = zeros( ComplexF64, tcd_dict[ (1, 1)])



  for left in 1:num
    for right in left+1:num

      tcd = tcd_dict[ (left, right)]
      occ .+= conj(phase_overlap[left]) .* tcd .* phase_overlap[right]
      occ .+= conj(phase_overlap[right]) .* conj.(tcd) .* phase_overlap[left]

    end 
  end 
  
  for idx in 1:num
    tcd = tcd_dict[ (idx, idx)]
    occ = occ.+ abs2( phase_overlap[idx]) * tcd
  end 


  return real.(occ), tcd_ex

end 


function get_eigen_tcd(tcd_dict, phase_overlap)

  num = minimum( map(maximum, collect(zip(collect(keys(tcd_dict))...))))

  #print(length(phase_overlap))
  #print(length(phase_overlap[1]))

  # we calculate the tcd between the time evolved GS state and every ex state of the system
  tcd_ex = [zeros( ComplexF64, length(tcd_dict[ (1, 1)])) for _ in 1:num]

  for left in 1:num
    for right in 1:num

      tcd = tcd_dict[ (left ,right)]
      tcd_ex[left] .+= conj(phase_overlap[right]) .* conj.(tcd)

    end 
  end 

  ## tcd_ex end ##

  return tcd_ex

end 






