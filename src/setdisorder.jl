"""Set the disorder in the system. Disorder has to be otherwise pre-generated, or set to 0"""
function setdisorder(disflag::Bool, L::Vector{Int})

  if disflag

    prefix = getworkdir()

    xraw = prefix * "disx"
    yraw = prefix * "disy"

    disx = readdlm(xraw)
    disy = readdlm(yraw)

  else
    
    disx = disy = zeros(1, prod(L))
  
  end 

  return disx, disy
end
