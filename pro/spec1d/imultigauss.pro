function imultigauss, logwave, pp, nline=nline, sigmares=sigmares
; jm08aug13nyu - adopted from existing code
    
    model = logwave*0.0
    for iline = 0L, nline-1L do begin
       p = pp[iline*3:iline*3+2]
       sigma = sqrt(sigmares[iline]^2.0 + p[2]^2.0)
       if (sigma ne 0.0) then begin
          term1 = exp(-0.5*(logwave - p[1])^2.0 / sigma^2.0)
          model = model + p[0]*term1/(sqrt(2.0*!pi)*sigma)
       endif
    endfor

return, model    
end    
