
function GS_wrapper()

    para = Dict(
        :L => 120,
        :N => 6,
        :sweepdim => 300,
        :sweepcnt => 50,
        :krylovdim => 8
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
        :sweepcnt => 50
        :QE => num,
        :QEen => energy,
        :ex => 25,
        :krylovdim => 8
    )

    QE(para)
end 


function QEdyna_wrapper()

    para = Dict(
        :L => 120,
        :N => 6,
        :sweepdim => 300,
        :sweepcnt => 50
        :QE => 2,
        :QEen => nothing,
        :dynamode => "left",
        :TEmethod => "TDVP",
        :TEdim => 150,
        :TEcutoff => 1E-8,
        :type => "Fermion",
        :screening => 0.0,
        :krylovdim => 8,
        :τ=>0.1,
        :QN=>true
    )

    QE_dynamic(para)
end 


function eigensolver_wrapper()

    additional_para = Dict{Any, Any}(
        :L => 12,
        :N => 6,
        :sweepdim => 300,
        :sweepcnt => 5,
        :krylovdim => 8,
    
      )

    eigensolver(additional_para)
end 