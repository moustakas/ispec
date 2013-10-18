pro istitch_2dspectra
; jm03dec21uofa

; cross-correlation parameters

    pmin = 5L
    pmax =  30L
    lags = lindgen((pmax-pmin)+1)+pmin
    nlags = n_elements(lags)
    npad = nlags/2L
    chi2array = fltarr(nlags)
    
; read the spectra
    
    fspec1 = 'fwcsra.2219.fits.gz'
    fspec2 = 'fwcsra.2220.fits.gz'

    spec1 = rd2dspec(fspec1) & spec2 = rd2dspec(fspec2)

; generate the large arrays that will be cross-correlated; pad by
; NLAGS/2 pixels; assume the images have the same dimensions
    
    npix = spec1.naxis2
    bignpix = 2*npix+2*npad
    bigindx = lindgen(bignpix)-npix
    bigrefspec = bigindx*0.0
    bigspec = bigindx*0.0

; positive pixels    

    refspec = im_normalize(total(spec1.image,1),/max)
    bigrefspec[npix+npad:2*npix+npad-1L] = refspec

    spec = im_normalize(total(spec2.image,1),/max)
    bigspec[npad:npad+npix-1L] = spec
    
    plot, bigindx, bigrefspec, ps=10, xsty=3, ysty=3
    oplot, bigindx, bigspec, ps=10

    for j = 0L, nlags-1L do chi2array[j] = total((shift(bigspec,lags[j])-bigrefspec)^2.0)
    plot, lags, chi2array, ps=10, xsty=3, ysty=3

stop    

return
end    
