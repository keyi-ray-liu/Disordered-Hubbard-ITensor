"""Set the disorder in the system. Disorder has to be otherwise pre-generated, or set to 0"""
function setdisorder(disflag::Bool, L)

  if disflag

    prefix = getworkdir()

    xraw = prefix * "disx"
    yraw = prefix * "disy"

    disx = readdlm(xraw)
    disy = readdlm(yraw)

  else
    
    if typeof(L) == Int

      disx = zeros((1, L))
      disy = zeros((1, L))

    else

      disx = zeros((1, L[1] * L[2]))
      disy = zeros((1, L[1] * L[2]))
    end 
  
  end 

  return disx, disy
end
