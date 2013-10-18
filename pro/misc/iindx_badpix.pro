function iindx_badpix, array, select_bit=select_bit, except_bit=except_bit, bit=bit
; jm03apr14uofa

    knownbits = fix([2,4,8,16,32])
    nbits = n_elements(knownbits)

    bit = fix(0)
    for i = 0L, nbits-1L do begin

       if total(array / knownbits[i] mod 2) gt 0.0 then bit = bit + knownbits[i]
    
    endfor

    indx = 0L
    
return, indx    
end
