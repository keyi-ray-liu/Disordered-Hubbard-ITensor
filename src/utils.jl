get_time(raw::String) = parse(Float64, SubString(raw, 1 + length(DYNA_STR), length(raw) - length(".h5")))

get_ex(raw::String, static_str::String) = parse(Int, SubString(raw, 1 + length(TEMP_tag) + length(static_str), length(raw) - length(".h5")))

get_start(raw::String) = parse(Int, SubString(raw, 1 + length("start"), length(raw) - length(".h5")))

get_dyna_files(; workflag = "") = sort( 
    filter(x-> !occursin("lasttime", x),
    filter(x->occursin(DYNA_STR,x), readdir(getworkdir(workflag))))
    , by=get_time)

shuffler(arr, random) = random ? shuffle(arr) : arr

get_static_files(static_str::String) = sort( filter(x->occursin( Regex(TEMP_tag * get_static_str("biasedchain") * ".*.h5"),x), readdir(getworkdir())), by= x -> get_ex(x, static_str ))

get_QE_ref_files() = sort( filter(x->(occursin(r"start.*h5",x) && !occursin(TEMP_tag, x)), readdir(getworkdir())), by=get_start)

positiveind(arr) = filter( x -> x!=0, [ arr[n] > 1e-5 ? n  : 0 for n in eachindex(arr)])

negativeind(arr) = filter( x -> x!=0, [ arr[n] < -1e-5 ? n  : 0 for n in eachindex(arr)])

FermionCondition(systype::String) = systype == "Fermion" ? 1 : 2


norm2(A,d=2) = sum(abs2,A,dims=d)

function FermionCondition(systype::String, t::Union{Number, Vector{T}}) where T <: Number

  if typeof(t) <: Number
    t = [ t for _ in 1:FermionCondition(systype)]
  end 
  @show t
  @assert length(t) == FermionCondition(systype) "t length does not match F/E systype"

  return t

end 


"""rewrite the function analogous to the original expect() function"""
function inner_product(L::MPS, R::MPS, opname::String; sites=1:length(L))

  @assert length(L) == length(R)
  s = siteinds(R)

  site_range = (sites isa AbstractRange) ? sites : collect(sites)
  Ns = length(site_range)


  ex = zeros(Complex, Ns)
  for (entry, j) in enumerate(site_range)
    
    Op =  op(opname, s[j])
    val = inner(L', apply(Op, R)) 
    ex[entry] = val

  end

  if sites isa Number
    return map(arr -> arr[1], ex)
  end
  return ex
end


function set_SD_parameters(s_coupling, d_coupling, contact_scaling, L::Int, config) 
  

  if config == "3x3"
      Lx = Ly = 3
      s_contacts = [ (s_coupling..., L + 1), (contact_scaling .* s_coupling..., L + 4), (s_coupling..., L + 7)]
      d_contacts = [ (d_coupling..., L + 3), (contact_scaling .* d_coupling..., L + 6), (d_coupling..., L + 9)]

  elseif config == "2x2"
      Lx = Ly = 2
      s_contacts = [ (s_coupling..., L + 1), (s_coupling..., L + 3)]
      d_contacts = [ (d_coupling..., L + 2), (d_coupling..., L + 4)]

  elseif config == "1x1"
      Lx = Ly = 1
      s_contacts = [(s_coupling..., L + 1)]
      d_contacts = [(d_coupling..., L + 1)]

  elseif config == "3x3-single" || config == "3x3-single-ring"
      Lx = Ly = 3
      s_contacts = [ (s_coupling..., L + 4)]
      d_contacts = [ (d_coupling..., L + 6)]

  elseif config == "3x3-19" 
      Lx = Ly = 3
      s_contacts = [ (s_coupling..., L + 1)]
      d_contacts = [ (d_coupling..., L + 9)]

  elseif config == "2x2-single"
      Lx = Ly = 2
      s_contacts = [ (s_coupling..., L + 1)]
      d_contacts = [ (d_coupling..., L + 4)]


  else
      error("Unrecognized config")
  end 
  
  return s_contacts, d_contacts, Lx, Ly
end 


function get_type_dict(systype)

    op_str = Dict(
      1 => systype == "Boson" ? "0" : "Emp",
      2 => systype == "Boson" ? "1" : systype == "Fermion" ? "Occ" : "Up",
      3 => systype == "Boson" ? "2" : "Dn",
      4 => systype == "Boson" ? "3" : "UpDn"
    )
  
    return op_str
  
end 

"""Set the work directory"""
function getworkdir( workflag )

  workdir = pwd() * "/work$(workflag)/"

  if !isdir(workdir)
    mkdir(workdir)

  end 

  return workdir


end 

vectomat( vec ) = mapreduce( permutedims, vcat, vec)


load_JSON(location) = JSON3.read(read(location, String))
#load_JSON(location) = JSON3.read(location,Union{Dict{String,Any},Array{Dict{String,Any}}})#JSON3.read(location, )

"""Wrapper function for the evaluation of the std of Hamiltonian"""
function variance(H::MPO, psi::MPS)
  @suppress begin 
    var = inner(H, psi, H, psi) - inner(psi', H, psi) ^ 2 
    return var
  end
end

check_ψ(output, workflag ) = isfile( getworkdir(workflag) * output * ".h5")


function load_ψ(output::String, workflag :: String; tag ="psi1")
  workdir = getworkdir(workflag )

  if length(output) > 3 && output[end-2:end] == ".h5"
    wf = h5open(workdir * output, "r")
  else
    wf = h5open( workdir * output * ".h5", "r")
  end 
  
  ψ = read(wf, tag, MPS)
  close(wf)

  return ψ
end 

function load_ψ(t::Float64, workflag :: String; tag ="psi1")
  workdir = getworkdir(workflag )
  wf = h5open( workdir * "tTDVP" * string(t) * ".h5", "r")
  ψ = read(wf, tag, MPS)
  close(wf)

  return ψ
end 

load_ψ(t::Int, workflag::String; tag ="psi1") = load_ψ( float(t), workflag; tag = tag)

function load_plasmon(output)
  ex = readdlm( getworkdir() * output * "ex")

  return ex[2] - ex[1]

end 

checkexist(output, workflag ) = isfile( getworkdir(workflag ) * output * ".h5")


# returns the precise site number of the ex QE site, beginning of the subsystem and end of subsystem, corresponding to site j
# function get_sys_loc(sys::QE_flat_SIAM, j::Int) 

#   if j <= left(sys)

#     mod = (j - 1) % (siteseach(sys) + QESITES)
#     qe_loc = j - mod 
#     chain_begin = qe_loc + QESITES
#     chain_end = chain_begin + siteseach(sys) - 1

#   elseif j > left(sys) + 1

#     mod = (j - 2) % (siteseach(sys) + QESITES)
#     chain_begin = j - mod 
#     chain_end = chain_begin + siteseach(sys) - 1
#     qe_loc = chain_end + 1


#   else
#     qe_loc = chain_begin = chain_end = -1
#   end 

#   return qe_loc, chain_begin, chain_end

# end 

# function gen_graph(sys::QE_G_SIAM)

#   g = NamedGraph()

#   for i in 1:legleft(sys) + legright(sys)
#     for j in 1:siteseach(sys) + QESITES
#       add_vertex!(g, (i, j))
#     end 
#   end 

#   center = (legleft(sys) + legright(sys) + 1, 1)
#   add_vertex!(g, center)

#   for i in 1:legright(sys) + legright(sys)

#     add_edge!(g, (center, (i, 1)))
#     for j in 1: siteseach(sys) +QESITES - 1
#       add_edge!(g, ((i, j), (i, j + 1)))
#     end 
#   end 

#   @show g
#   #@visualize g

#   return g
  

# end 

sitemap(::Systems, j::Int) = j
# sitemap(sys::Union{QE_G_SIAM, DPT_graph}, j) = sitemap(sys)[j]


"""maps the flattened index to graph index"""


# function get_sitemap(sys::QE_flat_SIAM)

#   d = Dict()

#   for j in 1:get_systotal(sys)

#     full = QESITES + siteseach(sys)

#     if j < left(sys) + 1
  
#       leg = div(j - 1, full) + 1
#       idx = full - (j - 1) % full
  
#     elseif j == left(sys) + 1
  
#       leg = legleft(sys) + legright(sys) + 1
#       idx = 1
  
#     else
  
#       j -= 1
#       leg = div(j - 1, full) + 1
#       idx = (j - 1) % full + 1
#       j += 1
  
#     end 
  
#     d[(leg, idx)] = j
#     d[j] = (leg, idx)


#   end 

#   @show d
#   return d
# end 


function gen_graph(sys::DPT_graph)

  g = NamedGraph()

  if typeof(sys.dpt) == DPT
    #L
    for j ∈ 1: L(sys) + R(sys)
      add_vertex!(g, (1, j))
    end 

    #dd
    for j ∈ 1:2
      add_vertex!(g, (2, j))
    end 


    for j ∈ 1:L(sys) + R(sys) - 1
      add_edge!(g, ((1, j), (1, j + 1)))
    end 

    add_edge!(g, ((2, 1), (1, L(sys))))
    add_edge!(g, ((2, 1), (2, 2)))

  else

    for j ∈ 1: L_contact(sys) - 1
      add_vertex!(g, (j, 1))
    end 

    for j ∈ 1 : 2 * couple_range(sys)
      add_vertex!(g, (L_contact(sys) , j))
    end 

    for j ∈ 1:2
      add_vertex!(g, (L_contact(sys) + 1, j))
    end 

    for j ∈ R_contact(sys) + 1:get_systotal(sys)
      
      leg = j - 2*couple_range(sys)
      add_vertex!(g, (leg, 1))
      add_edge!(g, ((leg, 1), (L_contact(sys), 2 * couple_range(sys))))
    end 

    for j ∈ 1: L_contact(sys) - 1
      add_edge!(g,((j, 1), (L_contact(sys), 1)))
    end 

    add_edge!(g, ((L_contact(sys) + 1, 1), (L_contact(sys) + 1, 2)))
    add_edge!(g, ((L_contact(sys) + 1, 1), (L_contact(sys) , couple_range(sys))))

    for j ∈ 1:2*couple_range(sys) - 1
      add_edge!(g, ((L_contact(sys), j), (L_contact(sys), j + 1)))
    end 

  end 

  #@show g
  #@visualize g

  return g
  

end 

function get_sitemap(sys::DPT)

  d = Dict()

  for j in 1:get_systotal(sys)

    if j <= L(sys)

      leg = 1
      idx = j

    elseif j <= L(sys) + 2
      leg = 2
      idx = j - L(sys)

    else
      leg = 1
      idx = j - 2

    end 
  
    d[(leg, idx)] = j
    d[j] = (leg, idx)


  end 

  @show d
  return d
end 


function get_sitemap(sys::DPT_mixed)

  d = Dict()

  for j in 1:get_systotal(sys)

    if j < L_contact(sys)

      leg = j
      idx = 1

    elseif j <= L(sys)
      leg = L_contact(sys)
      idx = j - L_contact(sys) + 1

    elseif j < R_begin(sys)
      leg = L_contact(sys) + 1
      idx = j - L(sys)

    elseif j <= R_contact(sys)
      leg = L_contact(sys)
      idx = j - L_contact(sys) - 1

    else
      leg = j - 2 * couple_range(sys) 
      idx = 1
    end 


  
    d[(leg, idx)] = j
    d[j] = (leg, idx)


  end 

  @show d
  return d
end 



function partial_contract(ψ::MPS, sites::Vector{Int})

  ψ = dense(ψ)
  s = siteinds(ψ)
  result = ITensor(1.)

  for j in eachindex(ψ)

    if j in sites
      result *= ψ[j]
    else 
      
      LHS = ITensor([1.0, 1.0],s[j])
      result *= LHS * ψ[j]

    end

  end

  #@show result
  return result
end 


function saveham(file, s, workflag)

  open( getworkdir(workflag) * "H" * file ,"w") do io
    write(io, string(s))
  end

end 


function gen_GS_scan()

  para_in = load_JSON( pwd() * "/qegaussian.json")
  full_size = get(para_in, "fullsize", 100)
  L = get(para_in, "L", 12)
  mode= get(para_in, "mode", "QE_two")

  # we need to define a singular site so that the later addition could proceed
  sites = siteinds("Fermion", full_size; conserve_qns=true)


  for start ∈ 1:L:(full_size - L + 1)
      para_in[ "chain_start" ] = start
      para_in[ "ex" ] = 1
      output = "start" * string(start)

      @info "GS scan start=$start"
      @show para_in

      solve_QE(; para_in = para_in, mode=mode, output=output, sites=sites)
  end 

end 


function get_tcd_gs()

  if !isfile( getworkdir() * "tcdgs.h5")

    if isempty(get_QE_ref_files())
        @warn "no gs raw files! generating GS scan wavefunctions"
        gen_GS_scan()
    end 

    @info "Summing reference state for TCD"

    @show get_QE_ref_files()
    gss = [ load_ψ(f) for f in get_QE_ref_files()]
    gs = add( gss..., maxdim = 128)

    normalize!(gs)

    h5open( getworkdir() * "tcdgs.h5", "w") do io
      write(io, "psi", gs)
    end 

  else

    @info "loading TCD gs"
    gs = load_ψ("tcdgs"; tag="psi")

  end 

  #@show expect(gs, "N")
  return gs

end 


function prev_res(workflag) 
  @info workflag

  if !isempty(Glob.glob( "$(LASTSTSTR).h5", getworkdir(workflag))) 

      @warn "prev state detected"

      try

        last_time = readdlm(getworkdir(workflag) * "tTDVPlasttime")[end]
        last_state = load_ψ( "$(LASTSTSTR)", workflag)
        
      catch 

        @warn "prev state broken, load back up state"
        last_time = readdlm(getworkdir(workflag) * "tTDVPbaktime")[end]
        last_state = load_ψ( "$(BACKSTR)", workflag)
        
      end

  else

      @warn "no prev state"
      last_time = -Inf
      last_state = MPS[]

  end 

  return last_time, last_state

end 



function ops_determiner(sys)

  if systype(sys) == "Fermion"
    ops = [ 
        ["Cdag", "C"],
        ["N", "N"]
    ]

  elseif typeof(sys) <: SD_array

      ops = [
          ["Cdagup", "Cup"],
          ["Cdagdn", "Cdn"],
          #["Nup", "Nup"],
          #["Ndn", "Ndn"]
      ]

  elseif typeof(sys) == DPT_avg
      ops = [
          ["Cdagup", "Cup"],
          #["Cdagdn", "Cdn"],
          #["Cdagup", "Cdn"],
          #["Cdagdn", "Cup"],
          ["Nup", "Nup"]
      ]
  end 

  return ops

end 


function get_time_control()

  timepara = load_JSON( pwd() * "/timepara.json")

  if haskey(timepara, "t1")
    τ1 = get(timepara, "t1", 0.25)
    τ2 = get(timepara, "t2", 2)
    fin1 = get(timepara, "fin1", 2.0)
    fin2 = get(timepara, "fin2", 100.0)

    timecontrol = TwoStage(τ1, τ2, fin1, fin2)
  else

    τ = get(timepara, "t", 0.25)
    fin = get(timepara, "fin", 100.0)

    timecontrol = OneStage(τ, fin)
  end 

  return timecontrol
end 



modifystate(ψ :: MPS, :: Identity, :: Systems) = ψ


function modifystate(ψ:: MPS, process :: LoadSource, sys :: Systems)

  Nop = systype(sys) == "Electron" ? "Ntot" : "N"
  leftinds, _ = reservoirmapping(sys)
  @show expect(ψ, Nop)

  s = siteinds(ψ)

  ops = systype(sys) == "Fermion" ? ["Cdag"] : ["Cdagup", "Cdagdn"]

  @show N = process.N

  sites_to_add = leftinds[  1 : N ]
  #reverse!(sites_to_add)
  for i in sites_to_add

    @show i, length(s), length(ψ)
    for opstring in ops

      operator = op(opstring, s[i])
      ψ = apply(operator, ψ)
    end 

  end 

  @show expect(ψ, Nop)
  return ψ

end 


function modifystate(ψ::MPS, process:: Perturbation, sys::Systems)

  ops = systype(sys) == "Fermion" ? ["Cdag"] : ["Cdagup", "Cdagdn"]
  s = siteinds(ψ)

  @show process

  for i in process.sites

    for opstring in ops
      operator = op(opstring, s[i])
      ψ = apply(operator, ψ)
    end 

  end 

  return ψ
end 



function get_mf(mode ;U = 0.0, mu = 0.0, L = 2, vs = 0.25 )
  
  if mode == "connectDD"
    #loop(; NLs = [L], vss= [vs], Us = [U], mus = [mu])
    val = []

    open(  pwd() * "/LOOP/U$(U)_vs$(vs)_mu$(mu)_L$(L)", "r") do IO
      val = readdlm(IO)
    end 

    @show val
    return val

  else
    return [0.0, 0.0]
  end 
end 





function two_site_rotate(ψ, α, ϕ)

  s = siteinds(ψ)
  M = zeros(ComplexF64, 4, 4)

  α = min(α, 1)
  M[2, 1] = sqrt(α)
  M[3, 1] = sqrt(1 - α) * exp( 1im * ϕ)
  #sqrt(0.7)
  #M[2, 2] = sqrt(0.3)

  a = op(M, s[end - 1], s[end])
  ψ = apply(a, ψ)

  return ψ
end 


function searchkey(s, key, func)

  s = split(s, "/")[end]
  strs = split(s, "_")

  for str in strs

    if str[1:length(key)] == key
      return func( str[ length(key) + 1 : end])
    end 

  end 

end 


function checkunique(fs, key, func)
  return sort([ searchkey(f, key, func) for f in fs])
end 

function fit()
    # @. linear(x, p) = p[1] * x + p[2]
    # p0 = [0.1, 0.5]

    # d1fit = curve_fit(linear, Us, d1s, p0)
    # @show coef(d1fit)

    # d1func(x) = linear(x, coef(d1fit))

end 