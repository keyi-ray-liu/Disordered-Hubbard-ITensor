
function GS_wrapper()

    para = Dict(
        :L => 120,
        :N => 6,
        :ex => 2,
        :sweepdim => 300,
        :sweepcnt => 50,
        :krylovdim => 8,
        :range => 10000
    )

    GSGap(para)
end 


function QE_wrapper(num, energy)

    num = parse(Int, num)
    energy = parse(Float64, energy)

    para = Dict(
        :L => 120,
        :N => 6,
        :sweepdim => 300,
        :sweepcnt => 50,
        :QE => num,
        :QEen => energy,
        :ex => 25,
        :krylovdim => 8
        :range => 10000
    )

    QE(para)
end 


function QEdyna_wrapper()

    para = Dict(
        :L => 120,
        :N => 6,
        :sweepdim => 300,
        :sweepcnt => 50,
        :QE => 2,
        :QEen => nothing,
        :dynamode => "left",
        :TEmethod => "TDVP",
        :TEdim => 150,
        :TEcutoff => 1E-8,
        :type => "Fermion",
        :screening => 0.0,
        :krylovdim => 8,
        :QN=>true
    )

    time_para = Dict(
        "τ" => 0.1,
        "start" => 0.1,
        "fin" => 100.0
    )

    QE_dynamic(para, time_para)
end 


function eigensolver_wrapper()

    GS_para = Dict{Any, Any}(
        :L => 120,
        :N => 60,
        :ex => 3,
        :sweepdim => 300,
        :sweepcnt => 25,
        :krylovdim => 8,
        :QN => true,
        :range => 10000
      )

    QE_internal_para = Dict(
        "QE_mul" => 1,
        "pl_level" => 2
    )

    QE_para = deepcopy(GS_para)

    QE_para[:ex] = 4
    QE_para[:N] = 60
    QE_para[:QE] = 2
    QE_para[:screening_qe] = 0.0

    dyna_para = deepcopy(QE_para)
    dyna_para[:dynamode]= "left"
    dyna_para[:TEmethod] = "eigen-occ"

    time_para = Dict(
        "τ" => 1.0,
        "start" => 1.0,
        "fin" => 1000.0
    )

    println(GS_para)
    println(QE_para)
    println(dyna_para)
    
    eigensolver(GS_para, QE_internal_para, QE_para, dyna_para, time_para)
end 