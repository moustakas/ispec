function read_standards, stdfile
; jm01jul12uofa
    
    nstds = n_elements(stdfile)
    spec = readfits(stdfile[0],spechead,/silent) ; read in the first spectrum
    specsize = (size(spec))[1]
    
    stdcube = fltarr(specsize,specsize,nstds)
    stdhead = ptr_arr(nstds)

    stdcube[*,*,0] = spec
    stdhead[0] = ptr_new(spechead)
    
    for k = 1L, nstds-1L do begin ; read in the remaining standards

       stdcube[*,*,k] = readfits(stdfile[k],spechead,/silent)
       stdhead[k] = ptr_new(spechead)       

    endfor

return, stdcube
end

pro fluxstar
; jm01jul11uofa

; stdfile is a list of filenames for all the standard stars that will
; be used to generate a sensitivity function
    
    stdfile = 'feige34.dat'
    path = '/home/ioannis/kennicutt/standards/'

    stdcube = read_standards(stdfile) ; read in all the standards into a spectrum cube

    extfile = 'kpnoextinct.dat'

    readcol, path+stdfile, stdwave, stdab, format='F,F', /silent   ; wavelength (A), ABmag
    readcol, path+extfile, extwave, extvals, format='F,F', /silent
 
; extracted standard star spectrum

    star = readfits('t.0131.0001.fits',stdhead,/silent) ; counts
    wave = make_wave(star,stdhead)                      ; wavelength

    stdname = sxpar(stdhead,'OBJECT')
    airmass = float(sxpar(stdhead,'AIRMASS'))
    exptime = float(sxpar(stdhead,'EXPTIME'))
    observatory, strn(sxpar(stdhead,'OBSERVAT')), obsinfo

    quadterp, extwave, extvals, wave, extinct ; interpolate the extinction wavelengths

    hcorrect = exp(-1.0*obsinfo.altitude/8300.0) ; height correction (what is the scale height of the troposphere!)
    factor = exp(extinct*hcorrect*airmass)       ; extinction correction

    quadterp, stdwave, stdab, wave, stdabcurve
    
    fstar = star * factor * 10D^(0.4*stdabcurve) / exptime ; erg/cm^2/s
    
stop    

return
end
