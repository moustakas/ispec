;+
; NAME:
;       PARSE_ISPECLINEFIT()
;
; PURPOSE:
;       Parse and write out the output from ISPECLINEFIT().
;
; INPUTS:
;       None.
;
; OPTIONAL INPUTS:
;       fitspecfile - output data structure from IFITSPEC
;       datapath    - data path to FITSPECFILE
;       outpath     - output path name (default DATAPATH)
;       root        - unique word to identify FITSPECFILE
;       ancillary   - ancillary data; to be used in place of PREPEND 
;       prepend     - data structure to prepend to the output from
;                     this routine (supersedes ancillary)
;       trimtags    - structure tags to trim from LINE (useful when
;                     passing PREPEND)
;       extra       - keywords for IUNRED_LINEDUST(),
;                     ICLASSIFICATION(), CONTINUUM_PROPERTIES(), and
;                     COMPUTE_LINELUMS() 
;
; KEYWORD PARAMETERS:
;       help        - print the syntax of this routine
;       silent      - do not print messages to STDOUT 
;       match       - match PREPEND and LINE before stitching 
;       write       - write the parsed data structure to OUTPATH 
;
; OUTPUTS:
;       line        - results of the emission-line fitting 
;
; OPTIONAL OUTPUTS:
;       linenodust  - same data structure as LINE but corrected for
;                     dust extinction using the Ha/Hb Balmer decrement 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 Mar 11, U of A - based on an earlier code 
;       jm04may19uofa - repair the [O III] 4959 line
;       jm04aug04uofa - repair the [N II] 6548 line; more error
;                       checking
;       jm04nov01uofa - store the Balmer emission-line fluxes and
;                       equivalent widths *uncorrected* for stellar
;                       absorption 
;       jm05jan06uofa - compute mass-to-light ratios separately from
;                       stellar masses
;       jm05jan07uofa - routine consolidated and speeded by
;                       suppressing extra calls to IUNRED_LINEDUST() 
;       jm05may17uofa - COMPUTE_STELLAR_MASS() is now done separately
;                       from this routine
;       jm05jul24uofa - remove TRIMTAGS from LINE rather than PREPEND 
;       jm05aug19uofa - added SYSERR optional input
;       jm05nov02uofa - added ANCILLARY optional input
;       jm06apr21uofa - added OUTPATH optional input
;       jm06nov18nyu  - added NOABUNDANCE keyword
;       jm08jan18nyu  - added SUFFIX optional input
;       jm08feb05nyu  - do not compute PROPERTIES or SFRS
;
; Copyright (C) 2004-2006, 2008, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function parse_ispeclinefit, fitspecfile=fitspecfile, datapath=datapath, $
  outpath=outpath, root=root, syserr=syserr, ancillary=ancillary, prepend=prepend1, $
  trimtags=trimtags, outfile=outfile, linenodust=linenodust, _extra=extra, $
  help=help, silent=silent, match=match, write=write, noabundance=noabundance, $
  suffix=suffix
    
    if keyword_set(help) then begin
       doc_library, 'parse_ispeclinefit'
       return, -1L
    endif 

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(outpath) eq 0L) then outpath = datapath
    if (n_elements(root) eq 0L) then root = ''
    if (n_elements(suffix) eq 0L) then suffix = ''
    if (n_elements(syserr) eq 0L) then syserr = 0.0
    if (syserr gt 0.0) then splog, 'WARNING: SYSERR > 0!'
    
    pushd, datapath

; restore the most recent line fit structure in datapath

    if n_elements(fitspecfile) eq 0L then begin
       allfiles = file_search(djs_filepath('*'+root+'*_specdata.fits.gz'),count=fcount)
       fitspecfile = allfiles[(reverse(sort(allfiles)))[0]]
    endif

    if file_test(fitspecfile) then begin
       if not keyword_set(silent) then splog, 'Reading '+datapath+fitspecfile+'.'
       line = mrdfits(fitspecfile,1,/silent)
    endif else begin
       if not keyword_set(silent) then splog, 'File '+datapath+fitspecfile+' not found.'
       return, -1L
    endelse

;   line = temporary(line[189163L:378325L])
;   line = temporary(line[378326L:567485L])
;   line = temporary(line[0L:189162L])

    if (n_elements(prepend1) ne 0L) then begin

       if keyword_set(match) then begin

; match in *both* structures simultaneously

          linegalaxy = strlowcase(strcompress(line.galaxy,/remove))
          pregalaxy = strlowcase(strcompress(prepend1.galaxy,/remove))

          match, linegalaxy, pregalaxy, lineindx, preindx
          line = line[lineindx]
          prepend = prepend1[preindx]

;         doit = match_string(linegalaxy,pregalaxy,index=index,/exact) ; <-- NOT GENERAL!!!
;         good = where(index ne -1L,ngood)
;         if (ngood ne 0L) then begin
;            prepend = prepend[index[good]]
;            line = line[good]
;         endif

       endif else begin ; assume the structures are already matched

          prepend = prepend1
          
       endelse

       if (n_elements(trimtags) ne 0L) then line = struct_trimtags(temporary(line),except=trimtags)
;      if (n_elements(trimtags) ne 0L) then prepend = struct_trimtags(prepend,except=trimtags)
       line = struct_addtags(prepend,temporary(line))
    
    endif
    nspec = n_elements(line)

    if (n_elements(prepend1) ne 0L) and (n_elements(ancillary) eq 0L) then ancillary = temporary(prepend)
    if (n_elements(ancillary) eq 0L) then ancillary = replicate({blank: ''},nspec)

    if (n_elements(ancillary) ne nspec) then begin
       splog, 'WARNING: ANCILLARY and LINE do not have the same number of elements!'
    endif

;   oldline = line

; add the systematic error in quadrature

    linename = strtrim(line[0].linename,2)
    nline = n_elements(linename)

    for iline = 0L, nline-1L do begin

; flux       
       
       true = tag_exist(line,linename[iline],index=fluxtag)
       flux = line.(fluxtag)[0,*]
       ferr = line.(fluxtag)[1,*]

       newferr = ferr

       pos = where(ferr gt 0.0,npos)
       if (npos ne 0L) then newferr[pos] = sqrt(ferr[pos]^2 + (flux[pos]*syserr/100.0)^2)

       line.(fluxtag) = [flux,newferr]
       
;;; ew - wrong!!
;;
;;       true = tag_exist(line,linename[iline]+'_ew',index=ewtag)
;;       ew = line.(ewtag)[0,*]
;;       ewerr = line.(ewtag)[1,*]
;;
;;       newewerr = ewerr
;;       
;;       pos = where(ewerr gt 0.0,npos)
;;       if (npos ne 0L) then newewerr[pos] = sqrt(ewerr[pos]^2 + (ew[pos]*syserr/100.0)^2)
;;
;;       line.(ewtag) = [ew,newewerr]
       
    endfor

    ew = 0 & ewerr = 0 & flux = 0 & ferr = 0 & newferr = 0 & newewerr = 0 & pos = 0 ; save memory
    
; repair the low S/N [O III] 4959 emission line by forcing it to be
; 1:3 with respect to [O III] 5007; append the original fluxes and
; errors to the output structure; if [O III] 5007 is an upper limit
; then use the same upper limit on [O III] 4959

    if tag_exist(line[0],'OIII_4959') and tag_exist(line[0],'OIII_5007') then begin

       oratio = 2.984 ; intrinsic ratio

       oiii = replicate({oiii_4959_orig: fltarr(2)},nspec)
       oiii.oiii_4959_orig = line.oiii_4959

       line.oiii_4959 = line.oiii_5007

       good = where(line.oiii_5007[1] gt 0.0,ngood)
       if (ngood ne 0L) then line[good].oiii_4959 = line[good].oiii_5007/oratio

       line = struct_addtags(temporary(line),temporary(oiii))
       
    endif
    
    if tag_exist(line[0],'OIII_4959_EW') and tag_exist(line[0],'OIII_5007_EW') then begin

       oratio = 2.984 ; intrinsic ratio

       oiii = replicate({oiii_4959_ew_orig: fltarr(2)},nspec)
       oiii.oiii_4959_ew_orig = line.oiii_4959_ew

       line.oiii_4959_ew = line.oiii_5007_ew

       good = where(line.oiii_5007_ew[1] gt 0.0,ngood)
       if (ngood ne 0L) then line[good].oiii_4959_ew = line[good].oiii_5007_ew/oratio

       line = struct_addtags(temporary(line),temporary(oiii))
       
    endif
    
    if tag_exist(line[0],'NII_6548') and tag_exist(line[0],'NII_6584') then begin

       nratio = 3.054

       nii = replicate({nii_6548_orig: fltarr(2)},nspec)
       nii.nii_6548_orig = line.nii_6548

       line.nii_6548 = line.nii_6584

       good = where(line.nii_6584[1] gt 0.0,ngood)
       if (ngood ne 0L) then line[good].nii_6548 = line[good].nii_6584/nratio

       upper = where(line.nii_6584[1] eq -3.0,nupper) ; upper limits
       if (nupper ne 0L) then line[upper].nii_6548[0] = line[upper].nii_6584[0]/nratio

       line = struct_addtags(temporary(line),temporary(nii))
       
    endif

    if tag_exist(line[0],'NII_6548_EW') and tag_exist(line[0],'NII_6584_EW') then begin

       nratio = 3.054

       nii = replicate({nii_6548_ew_orig: fltarr(2)},nspec)
       nii.nii_6548_ew_orig = line.nii_6548_ew

       line.nii_6548_ew = line.nii_6584_ew

       good = where(line.nii_6584_ew[1] gt 0.0,ngood)
       if (ngood ne 0L) then line[good].nii_6548_ew = line[good].nii_6584_ew/nratio

       line = struct_addtags(temporary(line),temporary(nii))
       
    endif

; add new tags giving the Balmer emission-line fluxes and equivalent
; widths *uncorrected* for stellar absorption

    if tag_exist(line[0],'BABS_H_ALPHA_EW') and tag_exist(line[0],'H_ALPHA') then begin
    
       halpha = replicate({h_alpha_uncor: fltarr(2), h_alpha_ew_uncor: fltarr(2)},nspec)

       halpha.h_alpha_uncor = line.h_alpha
       halpha.h_alpha_ew_uncor = line.h_alpha_ew

       good = where((line.h_alpha_ew[1] gt 0.0) and (line.babs_h_alpha_ew[1] gt 0.0),ngood)
       if (ngood ne 0L) then halpha[good].h_alpha_ew_uncor[0] = line[good].h_alpha_ew[0] - line[good].babs_h_alpha_ew[0]
       good = where((line.h_alpha[1] gt 0.0) and (line.babs_h_alpha[1] gt 0.0),ngood)
       if (ngood ne 0L) then halpha[good].h_alpha_uncor[0] = line[good].h_alpha[0] - line[good].babs_h_alpha[0]

;      plot, line[good].h_alpha_ew[0], halpha[good].h_alpha_ew_uncor[0]-line[good].h_alpha_ew[0], $
;        ps=4, /xlog, xsty=3, ysty=3
;      plot, line[good].h_alpha_ew[0], halpha[good].h_alpha_uncor[0]/line[good].h_alpha[0], $
;        ps=4, xsty=3, ysty=3, /xlog

       line = struct_addtags(temporary(line),temporary(halpha))
       
    endif
    
    if tag_exist(line[0],'BABS_H_BETA_EW') and tag_exist(line[0],'H_BETA') then begin
    
       hbeta = replicate({h_beta_uncor: fltarr(2), h_beta_ew_uncor: fltarr(2)},nspec)

       hbeta.h_beta_uncor = line.h_beta
       hbeta.h_beta_ew_uncor = line.h_beta_ew

       good = where((line.h_beta_ew[1] gt 0.0) and (line.babs_h_beta_ew[1] gt 0.0),ngood)
       if (ngood ne 0L) then hbeta[good].h_beta_ew_uncor[0] = line[good].h_beta_ew[0] - line[good].babs_h_beta_ew[0]
       good = where((line.h_beta[1] gt 0.0) and (line.babs_h_beta[1] gt 0.0),ngood)
       if (ngood ne 0L) then hbeta[good].h_beta_uncor[0] = line[good].h_beta[0] - line[good].babs_h_beta[0]

       line = struct_addtags(temporary(line),temporary(hbeta))
       
    endif
    
    if tag_exist(line[0],'BABS_H_GAMMA_EW') and tag_exist(line[0],'H_GAMMA') then begin
    
       hgamma = replicate({h_gamma_uncor: fltarr(2), h_gamma_ew_uncor: fltarr(2)},nspec)

       hgamma.h_gamma_uncor = line.h_gamma
       hgamma.h_gamma_ew_uncor = line.h_gamma_ew

       good = where((line.h_gamma_ew[1] gt 0.0) and (line.babs_h_gamma_ew[1] gt 0.0),ngood)
       if (ngood ne 0L) then hgamma[good].h_gamma_ew_uncor[0] = line[good].h_gamma_ew[0] - line[good].babs_h_gamma_ew[0]
       good = where((line.h_gamma[1] gt 0.0) and (line.babs_h_gamma[1] gt 0.0),ngood)
       if (ngood ne 0L) then hgamma[good].h_gamma_uncor[0] = line[good].h_gamma[0] - line[good].babs_h_gamma[0]

       line = struct_addtags(temporary(line),temporary(hgamma))
       
    endif
    
    if tag_exist(line[0],'BABS_H_DELTA_EW') and tag_exist(line[0],'H_DELTA') then begin
    
       hdelta = replicate({h_delta_uncor: fltarr(2), h_delta_ew_uncor: fltarr(2)},nspec)

       hdelta.h_delta_uncor = line.h_delta
       hdelta.h_delta_ew_uncor = line.h_delta_ew

       good = where((line.h_delta_ew[1] gt 0.0) and (line.babs_h_delta_ew[1] gt 0.0),ngood)
       if (ngood ne 0L) then hdelta[good].h_delta_ew_uncor[0] = line[good].h_delta_ew[0] - line[good].babs_h_delta_ew[0]
       good = where((line.h_delta[1] gt 0.0) and (line.babs_h_delta[1] gt 0.0),ngood)
       if (ngood ne 0L) then hdelta[good].h_delta_uncor[0] = line[good].h_delta[0] - line[good].babs_h_delta[0]

       line = struct_addtags(temporary(line),temporary(hdelta))
       
    endif

    good = 0 ; save memory
    
; compute continuum properties    

;   if not keyword_set(silent) then splog, 'Computing continuum properties.'
;   properties = continuum_properties(line,ancillary=ancillary,_extra=extra)
;   if (size(properties,/type) eq 8L) then line = struct_addtags(temporary(line),temporary(properties))

; apply classification criteria using the raw line fluxes 

    if not keyword_set(silent) then splog, 'Classifying.'
    iclass = iclassification(line,_extra=extra);,ratios=ratios)
    class = struct_trimtags(temporary(iclass),select=['*CLASS*'])
    if (size(class,/type) eq 8L) then line = struct_addtags(temporary(line),temporary(class))
    
; compute emission-line luminosities; (RELEGATED:) undefine ANCILLARY
; to save memory! 

    if not keyword_set(silent) then splog, 'Computing emission-line luminosities.'
    linelums = compute_linelums(line,ancillary=ancillary,_extra=extra)
;   linelums = compute_linelums(line,ancillary=temporary(ancillary),_extra=extra)
    if (size(linelums,/type) eq 8L) then line = struct_addtags(temporary(line),temporary(linelums))

; compute SFRS
    
;   if not keyword_set(silent) then splog, 'Computing emission-line SFRs.'
;   sfrs = compute_sfrs(line)
;   if (size(sfrs,/type) eq 8L) then line = struct_addtags(temporary(line),temporary(sfrs))

; compute strong-line abundances

    if (not keyword_set(noabundance)) then begin
       if not keyword_set(silent) then splog, 'Computing abundances.'
       abund = im_abundance(line,_extra=extra)
       if (size(abund,/type) eq 8L) then line = struct_addtags(temporary(line),temporary(abund))
    endif

; write out - save memory!
    
    if keyword_set(write) then begin

       if (n_elements(outfile) eq 0L) then outfile = $
         repstr(fitspecfile,'.fits.gz','_speclinefit'+suffix+'.fits')
       
       splog, 'Writing '+outpath+outfile+'.'
       mwrfits, line, outpath+outfile, /create
       spawn, ['gzip -f '+outpath+outfile], /sh
       
    endif

; get rid of tags that will be re-computed below    
    
    line = struct_trimtags(temporary(line),except=['SFR_*','ZSTRONG_*'])

; de-redden the line fluxes
    
    if not keyword_set(silent) then splog, 'De-reddening line fluxes.'
    linenodust = iunred_linedust(temporary(line),_extra=extra)
    if (size(linenodust,/type) ne 8L) then begin
       splog, 'De-reddening failed.'
       return, -1L
    endif

; compute SFRS 
    
;   if not keyword_set(silent) then splog, 'Computing emission-line SFRs.'
;   sfrsnodust = compute_sfrs(linenodust)
;   if (size(sfrsnodust,/type) eq 8L) then linenodust = struct_addtags(temporary(linenodust),temporary(sfrsnodust))
    
; compute strong-line abundances

    if (not keyword_set(noabundance)) then begin
       if not keyword_set(silent) then splog, 'Computing abundances.'
       abundnodust = im_abundance(linenodust,_extra=extra)
       if (size(abundnodust,/type) eq 8L) then linenodust = struct_addtags(temporary(linenodust),temporary(abundnodust))
    endif

; write out    
    
    if keyword_set(write) then begin

       outfile_nodust = repstr(outfile,'.fits','_nodust'+suffix+'.fits')
       
       splog, 'Writing '+outpath+outfile_nodust+'.'
       mwrfits, temporary(linenodust), outpath+outfile_nodust, /create
       spawn, ['gzip -f '+outpath+outfile_nodust], /sh
       
    endif

    popd
    
return, 1L
end    
