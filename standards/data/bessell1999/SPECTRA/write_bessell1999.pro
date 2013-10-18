pro write_bessell1999, debug=debug
; jm05jun01uofa
; write the Bessell (1999) standards (which are in mJy) into AB
; magnitudes

    rootpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='standards/data/bessell1999/SPECTRA')
    outpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='standards/data/bessell1999')
    ctiopath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='standards/obsolete/ctionewcal')

    pushd, rootpath

    flist = file_search('*.tab',count=fcount)
    outfile = 'bessell1999_'+repstr(flist,'.tab','.dat')
    star = strupcase(repstr(flist,'.tab',''))

    for i = 0L, fcount-1L do begin

       readcol, flist[i], wave, flux, skipline=7, format='F,F', /silent
       abmag = -2.5*alog10(flux*1D-29)-48.6

       openw, lun, outpath+outfile[i], /get_lun
       printf, lun, '# '+strtrim(star[i],2)
       for j = 0L, n_elements(wave)-1L do printf, lun, wave[j], abmag[j], 50.0, $
         format='(F7.1,3x,F6.3,3x,F4.1)'
       free_lun, lun

; to capture bugs, compare these magnitudes and wavelengths with
; CTIONEWCAL

       if keyword_set(debug) then begin
       
          readcol, ctiopath+repstr(outfile[i],'bessell1999','ctionewcal'), $
            cwave, cabmag, format='F,F', /silent

          if n_elements(wave) ne n_elements(cwave) then print, flist[i], 'Different number of wavelengths.'
;         print, 'Wavelength differences: ', flist[i], total(wave-cwave)
;         plot, wave, abmag, yrange=[max(abmag),min(abmag)], xsty=3, ysty=3, ps=4
;         oplot, cwave, cabmag, ps=7
          plot, wave, abmag-cabmag, yrange=[-0.6,0.1], $;xrange=[6500,10500], $
            xsty=3, ysty=3, ps=-4
          cc = get_kbrd(1)
          
       endif

    endfor

    popd
    
return
end
    
