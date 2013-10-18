;+
; NAME:
;       ISPEC_MAKELISTS
;
; PURPOSE:
;       Generate initial input file lists for ISPEC batch reduction.
;
; CALLING SEQUENCE:
;       ispec_makelists, root, outroot, datapath=, /overwrite, $
;          /all, /fitslist, /batchfile, /objlist, /stdlist, $
;          /tracelist, /seeinglist, /skylist, /crlist, /crsplits, $
;          /caliblist, /arclist, /gzext, _extra=extra
;
; INPUTS:
;       root    - unique string prefix (scalar or array) identifying a
;                 group of FITS files (see COMMENTS) [NNIGHTS]
;       outroot - unique output string suffix to append to the text
;                 file lists generated (see COMMENTS) [NNIGHTS]
;
; OPTIONAL INPUTS:
;       datapath  - I/O path
;       extra     - extra keywords for ISTARSEARCH() and
;                   WRITE_PARAMFILE
;       counter   - internal variable for WRITE_PARAMFILE
;
; KEYWORD PARAMETERS:
;       overwrite  - overwrite existing text file lists of the same 
;                    name without prompting
;       all        - generate all file lists (default)
;       fitslist   - only generate the FITS file list
;       batchfile  - only generate the IBATCH parameter files
;       objlist    - only generate the object file lists
;       stdlist    - only generate the standard star lists
;       tracelist  - only generate the distortion file lists
;       seeinglist - only generate the seeing file lists
;       skylist    - only generate the sky subtraction file lists and 
;                    the default SKYAPFILE
;       crsplits   - only generate the cosmic-ray split file lists  
;       crlist     - only generate the cosmic-ray file lists
;       caliblist  - only generate the calibration file lists
;       arclist    - only generate the arc lamp file lists
;       gzext      - write FITS.GZ extensions
;
; OUTPUTS:
;       See COMMENTS.
;
; OPTIONAL OUTPUTS:
;
; INTERNAL ROUTINES:
;       WRITE_PARAMFILE
;
; PROCEDURES USED:
;       IM_HMS2DEC(), MRDFITS(), STRUCT_APPEND, IM_DJS_ANGLE_MATCH(),
;       IFORAGE(), PRINT_STRUCT, PRECESS, CMSET_OP(), CWD(), SPLOG,
;       REMOVE, ISTARSEARCH()
;
; DATA FILES:
;       ${ISPEC_DIR}/standards/SUPPORTED
;
; TODO: 
;       Allow the user to customize the parameter file using
;       command-line optional inputs.
;
; COMMENTS:
;       ROOT should be like 'n1' or 'a.' to find all fits files with
;       file_search(root+'*.fits').  OUTROOT should be something like
;       '01dec20' to form, for example, the output file
;       'objlist_outroot.txt'.
;
;       This procedure will generate the following files:
;       'IBATCH_OUTROOT.TXT', 'OBJLIST_OUTROOT.TXT',
;       'STDLIST_OUTROOT.TXT', 'TRACELIST_OUTROOT.TXT',
;       'SKYLIST_OUTROOT.TXT', 'CRLIST_OUTROOT.TXT', and
;       'CALIBLIST_OUTROOT.TXT'. 
;
;       The standard stars are extracted from the object lists by
;       searching on position against the known database of
;       standards.
;
;       If the 'IBATCH_OUTROOT.TXT' parameter file exists and
;       ARCLIST=1 then also generate the WMAPNAME and ARCFILE files.  
;
;       Note that even if OVERWRITE is set, this routine will still
;       prompt before overwriting 'IBATCH_OUTROOT.TXT'.
;
;       This routine also generates a FITSLIST.TXT file which is
;       useful to print and have when reducing data. 
;
;       Arc lamps with no RA or DEC information in the header are not
;       written to ARCLIST and WMAPLIST.
;
; TODO:
;
; EXAMPLE:
;       Generate file lists and parameter files for a 3-night run that
;       took place in 2003 April 14-16.  The FITS files from each
;       night were written as 'a.10??.fits', 'a.11??.fits', and
;       'a.12??.fits'. 
;
;       IDL> ispec_makelists, ['a.10','a.11','a.12'], ['03apr14',$
;               '03apr15','03apr16'], datapath=datapath, /overwrite
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 October 6, U of A
;       jm03apr16uofa - re-written and documented, ISPEC v1.0.0
;       jm03apr22uofa - added additional error checking if no
;                       standards were observed on a particular night
;       jm03may27uofa - added CRSPLITS and ARCLIST keywords
;       jm03dec9uofa  - added GZIP, SEEING, and SKYLIST keywords;
;                       write default SKYAPFILE; wrote WRITE_PARAMFILE
;       jm04mar15uofa - syntactical updates to the file list comments 
;       jm04jun21uofa - bug fix: if only standard stars were observed
;                       then don't crash
;       jm04jun22uofa - do not generate SKYAPFILE lists
;       jm04nov05uofa - added SEARCHRAD optional input
;       jm05apr14uofa - replaced internal routine STARSEARCH() with
;                       the generalized ISTARSEARCH() and removed
;                       SEARCHRAD optional input; modifications to
;                       handle the new organization of the
;                       standard-star database; generalized
;                       WRITE_PARAMFILE 
;       jm05jun28uofa - prefixes changed to conform to v2.0 procedure 
;       arm05jul26uofa - bugs fixed: GZIP keyword replaced with GZEXT
;                        to avoid conflict with WRITE_PARAMFILE
;                        parameter, POPD added before return in
;                        recursive loop
;       jm07jun19nyu  - small bug fix (error checking) if no objects
;                       passed 
;
; Copyright (C) 2002-2005, John Moustakas
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

pro write_paramfile, paramfile, root=root, doplot=doplot, gzip=gzip, $
  biasfile=biasfile, domefile=domefile, domelist=domelist, skyflatfile=skyflatfile, $
  darkfile=darkfile, arcfile=arcfile, arclist=arclist, lampname=lampname, $
  flatname=flatname, illumname=illumname, wmapname=wmapname, tracename=tracename, $
  sensname=sensname, extfile=extfile, badpixfile=badpixfile, gain=gain, $
  rdnoise=rdnoise, trim=trim, overscan=overscan, pscale=pscale, $
  bsorder_flat=bsorder_flat, bsorder_sens=bsorder_sens, minwave_guess=minwave_guess, $
  minwave_out=minwave_out, dwave=dwave, sigclip=sigclip, objlim=objlim, sigfrac=sigfrac, $
  counter=counter

    if (n_elements(root) ne 0L) then if strcompress(root,/remove) ne '' then wroot = '_'+root

    name = [$
      'doplot          ',$
      'gzip            ',$
;     'saturation      ',$
      'biasfile        ',$
      'domefile        ',$
      'domelist        ',$
      'skyfile         ',$
      'darkfile        ',$
      'arcfile         ',$
      'arclist         ',$
      'lampname        ',$
      'flatname        ',$
      'illumname       ',$
      'wmapname        ',$
      'tracename       ',$
      'sensname        ',$
;     'stdpath         ',$
      'extfile         ',$
;     'skyspecfile     ',$
      'badpixfile      ',$
      'gain            ',$
      'rdnoise         ',$
      'trim            ',$
      'overscan        ',$
      'pscale          ',$
      'bsorder_flat    ',$
      'bsorder_sens    ',$
      'minwave_guess   ',$
      'minwave_out     ',$
      'dwave           ',$
      'sigclip         ',$
      'objlim          ',$
      'sigfrac         ']
    name = strupcase(name)
    nname = n_elements(name)

    value = [$
      ' 1                               ',$
      ' 0                               ',$
;     ' 0                               ',$
      'bias'+wroot+'.fits               ',$
      'dome'+wroot+'.fits               ',$
      ' 0                               ',$
      'skyflat'+wroot+'.fits            ',$
      '                                 ',$ ; darks not supported yet (jm05apr15uofa)
      'arclist'+wroot+'.txt             ',$
      ' 1                               ',$
      '                                 ',$
      'flat'+wroot+'.fits               ',$
      'illum'+wroot+'.fits              ',$
      'wmaplist'+wroot+'.txt            ',$
      'tracelist'+wroot+'.txt           ',$
      'sens'+wroot+'.fits               ',$
;     '                                 ',$
      '                                 ',$
;     '                                 ',$
      '                                 ',$
      '                                 ',$
      '                                 ',$
      '                                 ',$
      '                                 ',$
      '                                 ',$
      ' 50                              ',$
      ' 50                              ',$
      '                                 ',$
      '                                 ',$
      '                                 ',$
      ' 4.5                             ',$
      ' 2.0                             ',$
      ' 0.5                             ']
    nvalue = n_elements(value)

    comment = [$
      '# Boolean [1 means show plots]                                                 ',$
      '# compress [1 means gzip output FITS files]                                    ',$
;     '# pixel saturation value [electron] [0 means do not check for saturated pixels]',$
      '# master bias frame                                                            ',$
      '# master dome flat or file name                                                ',$
      '# Boolean [1 means DOMEFILE is a file name]                                    ',$
      '# sky flat                                                                     ',$
      '# master dark frame                                                            ',$
      '# arc lamp                                                                     ',$
      '# Boolean [1 means ARCFILE and WMAPNAME are file names]                        ',$
      '# type of comparison lamp                                                      ',$
      '# input/output flat field name                                                 ',$
      '# input/output illumation flat name                                            ',$
      '# input/output wavelength map name                                             ',$
      '# input/output trace structure                                                 ',$
      '# input/output sensitivity function name                                       ',$
;     '# subdirectory name of standard star files                                     ',$
      '# extinction file name [mag/airmass]                                           ',$
;     '# sky spectrum for wavelength cross-correlation                                ',$
      '# input bad pixel file                                                         ',$
      '# gain [electron/ADU]                                                          ',$
      '# read noise [electron]                                                        ',$
      '# trim region [zero-indexed]                                                   ',$
      '# overscan region [zero-indexed]                                               ',$
      '# spatial plate scale [arcsec/pixel]                                           ',$
      '# b-spline order to the flat fit                                               ',$
      '# b-spline order to the sensitivity function fit                               ',$
      '# approximate starting wavelength [Angstrom]                                   ',$
      '# output starting wavelength [Angstrom]                                        ',$
      '# approximate and output wavelength dispersion [Angstrom/pixel]                ',$
      '# cosmic ray clipping threshold [ILA_COSMIC]                                   ',$
      '# object clipping threshold [ILA_COSMIC]                                       ',$
      '# sigma fraction [ILA_COSMIC]                                                  ']
    ncomment = n_elements(comment)
    
    params = replicate({name: '', value: ' ', comment: ''},nname)
    params.name = name
    params.value = value
    params.comment = comment

    len = strlen(params[0].value)-1L

; required keywords
    
    if (n_elements(lampname) eq 0L) then splog, 'WARNING: required parameter LAMPNAME missing.' else $
      params[where(strlowcase(strcompress(name,/remove)) eq 'lampname')].value = $
        strtrim(lampname,2)+string(' ',format='(A'+string(len-strlen(strtrim(lampname,2)),format='(I0)')+')')
    if (n_elements(extfile) eq 0L) then splog, 'WARNING: required parameter EXTFILE missing.' else $
      params[where(strlowcase(strcompress(name,/remove)) eq 'extfile')].value = $
        strtrim(extfile,2)+string(' ',format='(A'+string(len-strlen(strtrim(extfile,2)),format='(I0)')+')')
    if (n_elements(gain) eq 0L) then splog, 'WARNING: required parameter GAIN missing.' else $
      params[where(strlowcase(strcompress(name,/remove)) eq 'gain')].value = $
        strtrim(gain,2)+string(' ',format='(A'+string(len-strlen(strtrim(gain,2)),format='(I0)')+')')
    if (n_elements(rdnoise) eq 0L) then splog, 'WARNING: required parameter RDNOISE missing.' else $
      params[where(strlowcase(strcompress(name,/remove)) eq 'rdnoise')].value = $
        strtrim(rdnoise,2)+string(' ',format='(A'+string(len-strlen(strtrim(rdnoise,2)),format='(I0)')+')')
    if (n_elements(trim) eq 0L) then splog, 'WARNING: required parameter TRIM missing.' else $
      params[where(strlowcase(strcompress(name,/remove)) eq 'trim')].value = $
        strtrim(trim,2)+string(' ',format='(A'+string(len-strlen(strtrim(trim,2)),format='(I0)')+')')
    if (n_elements(overscan) eq 0L) then splog, 'WARNING: required parameter OVERSCAN missing.' else $
      params[where(strlowcase(strcompress(name,/remove)) eq 'overscan')].value = $
        strtrim(overscan,2)+string(' ',format='(A'+string(len-strlen(strtrim(overscan,2)),format='(I0)')+')')
    if (n_elements(pscale) eq 0L) then splog, 'WARNING: required parameter PSCALE missing.' else $
      params[where(strlowcase(strcompress(name,/remove)) eq 'pscale')].value = $
        strtrim(pscale,2)+string(' ',format='(A'+string(len-strlen(strtrim(pscale,2)),format='(I0)')+')')
    if (n_elements(minwave_guess) eq 0L) then splog, 'WARNING: required parameter MINWAVE_GUESS missing.' else $
      params[where(strlowcase(strcompress(name,/remove)) eq 'minwave_guess')].value = $
        strtrim(minwave_guess,2)+string(' ',format='(A'+string(len-strlen(strtrim(minwave_guess,2)),format='(I0)')+')')
    if (n_elements(minwave_out) eq 0L) then splog, 'WARNING: required parameter MINWAVE_OUT missing.' else $
      params[where(strlowcase(strcompress(name,/remove)) eq 'minwave_out')].value = $
        strtrim(minwave_out,2)+string(' ',format='(A'+string(len-strlen(strtrim(minwave_out,2)),format='(I0)')+')')
    if (n_elements(dwave) eq 0L) then splog, 'WARNING: required parameter DWAVE missing.' else $
      params[where(strlowcase(strcompress(name,/remove)) eq 'dwave')].value = $
        strtrim(dwave,2)+string(' ',format='(A'+string(len-strlen(strtrim(dwave,2)),format='(I0)')+')')

; optional keywords assumed fixed for all parameter files

    if (n_elements(doplot) eq 1L) then params[where(strlowcase(strcompress(name,/remove)) eq 'doplot')].value = $
      strtrim(doplot,2)+string(' ',format='(A'+string(len-strlen(strtrim(doplot,2)),format='(I0)')+')')
    if (n_elements(gzip) eq 1L) then params[where(strlowcase(strcompress(name,/remove)) eq 'gzip')].value = $
      strtrim(gzip,2)+string(' ',format='(A'+string(len-strlen(strtrim(gzip,2)),format='(I0)')+')')
    if (n_elements(domelist) eq 1L) then params[where(strlowcase(strcompress(name,/remove)) eq 'domelist')].value = $
      strtrim(domelist,2)+string(' ',format='(A'+string(len-strlen(strtrim(domelist,2)),format='(I0)')+')')
    if (n_elements(arclist) eq 1L) then params[where(strlowcase(strcompress(name,/remove)) eq 'arclist')].value = $
      strtrim(arclist,2)+string(' ',format='(A'+string(len-strlen(strtrim(arclist,2)),format='(I0)')+')')
    if (n_elements(badpixfile) eq 1L) then params[where(strlowcase(strcompress(name,/remove)) eq 'badpixfile')].value = $
      strtrim(badpixfile,2)+string(' ',format='(A'+string(len-strlen(strtrim(badpixfile,2)),format='(I0)')+')')
    if (n_elements(bsorder_flat) eq 1L) then params[where(strlowcase(strcompress(name,/remove)) eq 'bsorder_flat')].value = $
      strtrim(bsorder_flat,2)+string(' ',format='(A'+string(len-strlen(strtrim(bsorder_flat,2)),format='(I0)')+')')
    if (n_elements(bsorder_sens) eq 1L) then params[where(strlowcase(strcompress(name,/remove)) eq 'bsorder_sens')].value = $
      strtrim(bsorder_sens,2)+string(' ',format='(A'+string(len-strlen(strtrim(bsorder_sens,2)),format='(I0)')+')')
    if (n_elements(sigclip) eq 1L) then params[where(strlowcase(strcompress(name,/remove)) eq 'sigclip')].value = $
      strtrim(sigclip,2)+string(' ',format='(A'+string(len-strlen(strtrim(sigclip,2)),format='(I0)')+')')
    if (n_elements(objlim) eq 1L) then params[where(strlowcase(strcompress(name,/remove)) eq 'objlim')].value = $
      strtrim(objlim,2)+string(' ',format='(A'+string(len-strlen(strtrim(objlim,2)),format='(I0)')+')')
    if (n_elements(sigfrac) eq 1L) then params[where(strlowcase(strcompress(name,/remove)) eq 'sigfrac')].value = $
      strtrim(sigfrac,2)+string(' ',format='(A'+string(len-strlen(strtrim(sigfrac,2)),format='(I0)')+')')

; optional keywords potentially different for each parameter file

    if (n_elements(biasfile) ge 1L) then begin
       if (n_elements(biasfile) eq 1L) then indx = 0L else indx = counter
       if (indx gt n_elements(biasfile)-1L) then splog, 'BIASFILE has the wrong number of elements.' else $
         params[where(strlowcase(strcompress(name,/remove)) eq 'biasfile')].value = $
           strtrim(biasfile[indx],2)+string(' ',format='(A'+string(len-strlen(strtrim(biasfile[indx],2)),format='(I0)')+')')
    endif

    if (n_elements(domefile) ge 1L) then begin
       if (n_elements(domefile) eq 1L) then indx = 0L else indx = counter
       if (indx gt n_elements(domefile)-1L) then splog, 'DOMEFILE has the wrong number of elements.' else $
         params[where(strlowcase(strcompress(name,/remove)) eq 'domefile')].value = $
           strtrim(domefile[indx],2)+string(' ',format='(A'+string(len-strlen(strtrim(domefile[indx],2)),format='(I0)')+')')
    endif

    if (n_elements(skyflatfile) ge 1L) then begin
       if (n_elements(skyflatfile) eq 1L) then indx = 0L else indx = counter
       if (indx gt n_elements(skyflatfile)-1L) then splog, 'SKYFLATFILE has the wrong number of elements.' else $
         params[where(strlowcase(strcompress(name,/remove)) eq 'skyfile')].value = $
           strtrim(skyflatfile[indx],2)+string(' ',format='(A'+string(len-strlen(strtrim(skyflatfile[indx],2)),format='(I0)')+')')
    endif

    if (n_elements(darkfile) ge 1L) then begin
       if (n_elements(darkfile) eq 1L) then indx = 0L else indx = counter
       if (indx gt n_elements(darkfile)-1L) then splog, 'DARKFILE has the wrong number of elements.' else $
         params[where(strlowcase(strcompress(name,/remove)) eq 'darkfile')].value = $
           strtrim(darkfile[indx],2)+string(' ',format='(A'+string(len-strlen(strtrim(darkfile[indx],2)),format='(I0)')+')')
    endif

    if (n_elements(arcfile) ge 1L) then begin
       if (n_elements(arcfile) eq 1L) then indx = 0L else indx = counter
       if (indx gt n_elements(arcfile)-1L) then splog, 'ARCFILE has the wrong number of elements.' else $
         params[where(strlowcase(strcompress(name,/remove)) eq 'arcfile')].value = $
           strtrim(arcfile[indx],2)+string(' ',format='(A'+string(len-strlen(strtrim(arcfile[indx],2)),format='(I0)')+')')
    endif

    if (n_elements(flatname) ge 1L) then begin
       if (n_elements(flatname) eq 1L) then indx = 0L else indx = counter
       if (indx gt n_elements(flatname)-1L) then splog, 'FLATNAME has the wrong number of elements.' else $
         params[where(strlowcase(strcompress(name,/remove)) eq 'flatname')].value = $
           strtrim(flatname[indx],2)+string(' ',format='(A'+string(len-strlen(strtrim(flatname[indx],2)),format='(I0)')+')')
    endif

    if (n_elements(illumname) ge 1L) then begin
       if (n_elements(illumname) eq 1L) then indx = 0L else indx = counter
       if (indx gt n_elements(illumname)-1L) then splog, 'ILLUMNAME has the wrong number of elements.' else $
         params[where(strlowcase(strcompress(name,/remove)) eq 'illumname')].value = $
           strtrim(illumname[indx],2)+string(' ',format='(A'+string(len-strlen(strtrim(illumname[indx],2)),format='(I0)')+')')
    endif

    if (n_elements(wmapname) ge 1L) then begin
       if (n_elements(wmapname) eq 1L) then indx = 0L else indx = counter
       if (indx gt n_elements(wmapname)-1L) then splog, 'WMAPNAME has the wrong number of elements.' else $
         params[where(strlowcase(strcompress(name,/remove)) eq 'wmapname')].value = $
           strtrim(wmapname[indx],2)+string(' ',format='(A'+string(len-strlen(strtrim(wmapname[indx],2)),format='(I0)')+')')
    endif

    if (n_elements(tracename) ge 1L) then begin
       if (n_elements(tracename) eq 1L) then indx = 0L else indx = counter
       if (indx gt n_elements(tracename)-1L) then splog, 'TRACENAME has the wrong number of elements.' else $
         params[where(strlowcase(strcompress(name,/remove)) eq 'tracename')].value = $
           strtrim(tracename[indx],2)+string(' ',format='(A'+string(len-strlen(strtrim(tracename[indx],2)),format='(I0)')+')')
    endif

    if (n_elements(sensname) ge 1L) then begin
       if (n_elements(sensname) eq 1L) then indx = 0L else indx = counter
       if (indx gt n_elements(sensname)-1L) then splog, 'SENSNAME has the wrong number of elements.' else $
         params[where(strlowcase(strcompress(name,/remove)) eq 'sensname')].value = $
           strtrim(sensname[indx],2)+string(' ',format='(A'+string(len-strlen(strtrim(sensname[indx],2)),format='(I0)')+')')
    endif

; write out
    
    openw, lun, paramfile, /get_lun
    struct_print, params, /no_head, lun=lun
    free_lun, lun

return
end

pro ispec_makelists, root, outroot, datapath=datapath, overwrite=overwrite, all=all, $
  fitslist=fitslist, batchfile=batchfile, objlist=objlist, stdlist=stdlist, $
  tracelist=tracelist, seeinglist=seeinglist, skylist=skylist, crlist=crlist, $
  crsplits=crsplits, caliblist=caliblist, arclist=arclist, gzext=gzext, _extra=extra
    
    if n_params() ne 2L then begin
       print, 'Syntax - ispec_makelists, root, outroot, datapath=, /overwrite, $'
       print, '   /all, /fitslist, /batchfile, /objlist, /stdlist, /tracelist, $'
       print, '   /seeinglist, /skylist, /crlist, /crsplits, /caliblist, /arclist, $'
       print, '   /gzext, _extra=extra'
       return
    endif
    
    nroot = n_elements(root)
    noutroot = n_elements(outroot)

    if nroot ne noutroot then begin
       print, 'ROOT and OUTROOT must have the same number of elements.'
       return
    endif

    if n_elements(datapath) eq 0L then datapath = cwd()
    etcpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    
    pushd, datapath

; vectorize
    
    if nroot gt 1L then begin

       for k = 0L, nroot-1L do ispec_makelists, root[k], outroot[k], datapath=datapath, $
         overwrite=overwrite, all=all, fitslist=fitslist, batchfile=batchfile, $
         objlist=objlist, stdlist=stdlist, tracelist=tracelist, seeinglist=seeinglist, $
         skylist=skylist, crlist=crlist, crsplits=crsplits, caliblist=caliblist, $
         arclist=arclist, gzext=gzext, counter=k, _extra=extra
       popd
       return
       
    endif

; make the default: generate all file lists    
    
    if (n_elements(all) eq 0L) and (n_elements(fitslist) eq 0L) and $
      (n_elements(batchfile) eq 0L) and (n_elements(objlist) eq 0L) and $
      (n_elements(stdlist) eq 0L) and (n_elements(tracelist) eq 0L) and $
      (n_elements(seeinglist) eq 0L) and (n_elements(skylist) eq 0L) and $
      (n_elements(crlist) eq 0L) and (n_elements(crsplits) eq 0L) and $
      (n_elements(caliblist) eq 0L) and (n_elements(arclist) eq 0L) then $
      all = 1L
    
; gather useful header information on each spectrum and segregate
; "objects" from "other" calibration data

    flist = file_search(root+'*.fits*',count=fcount)
    info = iforage(flist)
    objects = where(strmatch(info.type,'*object*',/fold) eq 1B,nobject,comp=other,ncomp=nother)

    if nobject ne 0L then objinfo = info[objects]
    if nother ne 0L then otherinfo = info[other]
    
; write the FITS file list for all objects

    if keyword_set(all) or keyword_set(fitslist) then begin

       write = 0L

       fitslistfile = outroot+'_fitslist.txt'
       if file_test(fitslistfile,/regular) and (not keyword_set(overwrite)) then begin
          print, 'Overwrite existing file '+fitslistfile+' [Y/N]? ', format='(A,$ )'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then write = 1L
          print
       endif else write = 1L

;      objmax = strmax(strlen(info.object))

       if write then begin
          openw, lun1, fitslistfile, /get_lun
          print_struct, info, ['FILE','TYPE','OBJECT'], title, outstr, /strings
          printf, lun1, '# '+title
          for k = 0L, n_elements(outstr[0,*])-1L do printf, lun1, outstr[*,k], format='(3A)'
          free_lun, lun1
       endif
    endif
       
; find the standard stars 

    if (nobject ne 0L) then begin

       starinfoindx = istarsearch(objinfo.ra,objinfo.dec,epoch=objinfo.epoch,$
         objindx=starindx,_extra=extra)
       nonstarindx = -1L

       if (starindx[0] ne -1L) then begin

          nonstarindx = lindgen(nobject)
          if (nobject ne n_elements(starindx)) then $ ; jm04jun21uofa
            remove, starindx, nonstarindx else $ 
              nonstarindx = -1L

          sinfo = objinfo[starindx]
          sinfo.file = 'sr'+strcompress(sinfo.file,/remove)
          if keyword_set(gzext) then sinfo.file = sinfo.file+'.gz'
       endif else splog, 'No standard stars observed '+root+', '+outroot+'.'

    endif else starindx = -1L

; initialize the IBATCH parameter file

    if keyword_set(all) or keyword_set(batchfile) then begin
    
       write = 0L
       
       paramfile = 'ibatch_'+outroot+'.txt'
       if file_test(paramfile,/regular) then begin
;      if file_test(paramfile,/regular) and (not keyword_set(overwrite)) then begin
          print, 'Overwrite existing file '+paramfile+' [Y/N]? ', format='(A,$ )'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then write = 1L
          print
       endif else write = 1L

       if write then write_paramfile, paramfile, root=outroot, counter=counter, _extra=extra
;      if write then spawn, ['/bin/cp -f '+etcpath+'ibatch_template.txt ./'+paramfile], /sh

    endif
       
; write the OBJLIST file

    if keyword_set(all) or keyword_set(objlist) then begin

       write = 0L

       objfile = 'objlist_'+outroot+'.txt'
       if file_test(objfile,/regular) and (not keyword_set(overwrite)) then begin
          print, 'Overwrite existing file '+objfile+' [Y/N]? ', format='(A,$ )'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then write = 1L
          print
       endif else write = 1L
       
       if write then begin
          openw, lun1, objfile, /get_lun
          printf, lun1, '# File list for ICCDPROC.  Do not include bias frames, dome flats'
          printf, lun1, '# or darks but do include arc lamps and, optionally, sky flats.'
          print_struct, info, ['FILE','OBJECT','TYPE'], title, outstr, /strings
          printf, lun1, '# '+title
          for k = 0L, n_elements(outstr[0,*])-1L do printf, lun1, '  ', outstr[*,k], format='(A2,3A)'
          free_lun, lun1
       endif

    endif
       
; write the STDLIST file (standard star list)

    if (keyword_set(all) or keyword_set(stdlist)) and (starindx[0] ne -1L) then begin

       write = 0L
       
       stdfile = 'stdlist_'+outroot+'.txt'
       if file_test(stdfile,/regular) and (not keyword_set(overwrite)) then begin
          print, 'Overwrite existing file '+stdfile+' [Y/N]? ', format='(A,$ )'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then write = 1L
          print
       endif else write = 1L
       
       if write then begin
          openw, lun1, stdfile, /get_lun
          printf, lun1, '# Standard star file list. '
          print_struct, sinfo, ['FILE','OBJECT','TYPE','AIRMASS','APERTURE'], title, outstr, /strings
;         print_struct, sinfo, ['FILE','OBJECT','TYPE','AIRMASS','APERTURE','PA','PARANGLE'], title, outstr, /strings
          printf, lun1, '# '+title
          for k = 0L, n_elements(outstr[0,*])-1L do printf, lun1, '  ', outstr[*,k], format='(A2,7A)'
          free_lun, lun1
       endif

    endif
       
; write the TRACELIST file (use the standard star list)

    if (keyword_set(all) or keyword_set(tracelist)) and (starindx[0] ne -1L) then begin

       write = 0L
       tinfo = objinfo[starindx]
       tinfo.file = 'r'+strcompress(tinfo.file,/remove)
       if keyword_set(gzext) then tinfo.file = tinfo.file+'.gz'
       
       tracefile = 'tracelist_'+outroot+'.txt'
       if file_test(tracefile,/regular) and (not keyword_set(overwrite)) then begin
          print, 'Overwrite existing file '+tracefile+' [Y/N]? ', format='(A,$ )'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then write = 1L
          print
       endif else write = 1L
       
       if write then begin
          openw, lun1, tracefile, /get_lun
          printf, lun1, '# Spatial distortion file list. '
          print_struct, tinfo, ['FILE','OBJECT','TYPE','AIRMASS','APERTURE'], title, outstr, /strings
          printf, lun1, '# '+title
          for k = 0L, n_elements(outstr[0,*])-1L do printf, lun1, '  ', outstr[*,k], format='(A2,5A)'
          free_lun, lun1
       endif

    endif
       
; write the SEEINGLIST file (use the standard star list)

    if (keyword_set(all) or keyword_set(seeinglist)) and (starindx[0] ne -1L) then begin

       write = 0L
       sinfo.file = 'w'+strcompress(sinfo.file,/remove)
       
       seeingfile = 'seeing_'+outroot+'.txt'
       if file_test(seeingfile,/regular) and (not keyword_set(overwrite)) then begin
          print, 'Overwrite existing file '+seeingfile+' [Y/N]? ', format='(A,$ )'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then write = 1L
          print
       endif else write = 1L
       
       if write then begin
          openw, lun1, seeingfile, /get_lun
          printf, lun1, '# Seeing file list. '
          print_struct, sinfo, ['FILE','OBJECT','TYPE','AIRMASS','APERTURE'], title, outstr, /strings
;         print_struct, sinfo, ['FILE','OBJECT','TYPE','AIRMASS','APERTURE','PA','PARANGLE'], title, outstr, /strings
          printf, lun1, '# '+title
          for k = 0L, n_elements(outstr[0,*])-1L do printf, lun1, '  ', outstr[*,k], format='(A2,7A)'
          free_lun, lun1
       endif

    endif
       
; write the SKYLIST and the SKYAPFILE files

    if keyword_set(all) or keyword_set(skylist) then begin

       write = 0L
       skinfo = objinfo
       skinfo.file = 'cr'+strcompress(skinfo.file,/remove)
       if keyword_set(gzext) then skinfo.file = skinfo.file+'.gz'

; group the standard stars at the top of these files
       
       if (starindx[0] ne -1L) and (nonstarindx[0] ne -1L) then begin ; jm04jun21uofa
          skinfo = struct_append(skinfo[starindx],skinfo[nonstarindx])
       endif

       skyfile = 'skylist_'+outroot+'.txt'
       if file_test(skyfile,/regular) and (not keyword_set(overwrite)) then begin
          print, 'Overwrite existing file '+skyfile+' [Y/N]? ', format='(A,$ )'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then write = 1L
          print
       endif else write = 1L
       
       if write then begin
          openw, lun1, skyfile, /get_lun
          printf, lun1, '# Sky subtraction (ISKYSUBTRACT2D) file list. Do not include biases,'
          printf, lun1, '# dome flats, darks, arc lamps, and (optionally) sky flats.'
          print_struct, skinfo, ['FILE','OBJECT','TYPE'], title, outstr, /strings
          printf, lun1, '# '+title
          for k = 0L, n_elements(outstr[0,*])-1L do printf, lun1, '  ', outstr[*,k], format='(A2,3A)'
          free_lun, lun1
       endif

;      skinfo.object = strcompress(skinfo.object,/remove)
;      skyapfile = 'skyaplist_'+outroot+'.txt'
;      if file_test(skyapfile,/regular) and (not keyword_set(overwrite)) then begin
;         print, 'Overwrite existing file '+skyapfile+' [Y/N]? ', format='(A,$ )'
;         cc = get_kbrd(1)
;         if strupcase(cc) eq 'Y' then write = 1L
;         print
;      endif else write = 1L
;      
;      if write then begin
;         temp = replicate({skylower: 0.0, skyupper: 0.0, skyaperturelo: 20.0, $
;           skyaperturehi: 20.0, norder_sky: 2},nobject)
;         skinfo = struct_addtags(skinfo,temp)
;         
;         openw, lun1, skyapfile, /get_lun
;         printf, lun1, '# Sky apertures (ISKYSUBTRACT2D) file list.  This file should mirror '+skyfile+'.'
;         print_struct, skinfo, ['FILE','OBJECT','SKYLOWER','SKYUPPER','SKYAPERTURELO','SKYAPERTUREHI','NORDER_SKY'], $
;           title, outstr, /strings
;         printf, lun1, '# '+title
;         for k = 0L, n_elements(outstr[0,*])-1L do printf, lun1, '  ', outstr[*,k], format='(A2,7A)'
;         free_lun, lun1
;      endif

   endif
       
; write the CRLIST file

    if keyword_set(all) or keyword_set(crlist) then begin

       write = 0L
       cinfo = objinfo
       cinfo.file = 'r'+strcompress(cinfo.file,/remove)
       if keyword_set(gzext) then cinfo.file = cinfo.file+'.gz'
       
       crfile = 'crlist_'+outroot+'.txt'
       if file_test(crfile,/regular) and (not keyword_set(overwrite)) then begin
          print, 'Overwrite existing file '+crfile+' [Y/N]? ', format='(A,$ )'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then write = 1L
          print
       endif else write = 1L
       
       if write then begin
          openw, lun1, crfile, /get_lun
          printf, lun1, '# Cosmic ray rejections (LA_COSMIC) file list. Do not include biases,'
          printf, lun1, '# dome flats, darks, sky flats, arc lamps, or (optionally) standard stars.' 
          print_struct, cinfo, ['FILE','OBJECT','TYPE'], title, outstr, /strings
          printf, lun1, '# '+title
          for k = 0L, n_elements(outstr[0,*])-1L do printf, lun1, '  ', outstr[*,k], format='(A2,3A)'
          free_lun, lun1
       endif

    endif
       
; write the CRSPLITS file

    if keyword_set(all) or keyword_set(crsplits) then begin

       write = 0L

       cinfo = objinfo
       cinfo.file = 'r'+strcompress(cinfo.file,/remove)
       if keyword_set(gzext) then cinfo.file = cinfo.file+'.gz'
       
       crfile = 'crsplits_'+outroot+'.txt'
       if file_test(crfile,/regular) and (not keyword_set(overwrite)) then begin
          print, 'Overwrite existing file '+crfile+' [Y/N]? ', format='(A,$ )'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then write = 1L
          print
       endif else write = 1L

       if write then begin

          openw, lun1, crfile, /get_lun
          printf, lun1, '# Cosmic ray splits (ICRCOMBINE) file list.  Verify against observations logs. '
          printf, lun1, '# '

          bigindx = lindgen(nobject)
          for k = 0L, nobject-1L do begin
             
             doit = match_string(cinfo[k].object,cinfo[bigindx].object,index=match)
             nmatch = n_elements(match)
             if (match[0] ne -1L) and (nmatch gt 1L) then begin

                len = strn(max(strlen(cinfo[bigindx[match]].file)))
;               print, cinfo[bigindx[match]].file, format='('+strjoin(replicate('A'+len,nmatch),',1x,')+')'
                printf, lun1, cinfo[bigindx[match]].file, format='('+strjoin(replicate('A'+len,nmatch),',1x,')+')'
                if (n_elements(bigindx) ne nmatch) then remove, match, bigindx
                
             endif
             
          endfor

          free_lun, lun1

       endif

    endif 
       
; write the CALIBLIST file

    if keyword_set(all) or keyword_set(caliblist) then begin

       write = 0L
       cinfo = objinfo
       cinfo.file = 'scr'+strcompress(cinfo.file,/remove)
       if keyword_set(gzext) then cinfo.file = cinfo.file+'.gz'
       
       calibfile = 'caliblist_'+outroot+'.txt'
       if file_test(calibfile,/regular) and (not keyword_set(overwrite)) then begin
          print, 'Overwrite existing file '+calibfile+' [Y/N]? ', format='(A,$ )'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then write = 1L
          print
       endif else write = 1L
       
       if write then begin
          openw, lun1, calibfile, /get_lun
          printf, lun1, '# Calibration file list.  '+$
            'Do not include biases, dome flats, darks, sky flats, or arc lamps.'
          print_struct, cinfo, ['FILE','OBJECT','TYPE'], title, outstr, /strings
          printf, lun1, '# '+title
          for k = 0L, n_elements(outstr[0,*])-1L do printf, lun1, '  ', outstr[*,k], format='(A2,3A)'
          free_lun, lun1
       endif

    endif
       
; write the ARCLIST file

    if keyword_set(all) or keyword_set(arclist) then begin

       write = 0L

       arcs = where((strmatch(otherinfo.type,'*comp*',/fold) eq 1B) and $
         (otherinfo.ra ne '') and (otherinfo.dec ne ''),narcs)
       if narcs eq 0L then splog, 'ARCLIST: No arc lamps in file list.' else begin
       
          ainfo = otherinfo[arcs]
          ainfo.file = 'r'+strcompress(ainfo.file,/remove)
          if keyword_set(gzext) then ainfo.file = ainfo.file+'.gz'

          arcfile = 'arclist_'+outroot+'.txt'
          if file_test(arcfile,/regular) and (not keyword_set(overwrite)) then begin
             print, 'Overwrite existing file '+arcfile+' [Y/N]? ', format='(A,$ )'
             cc = get_kbrd(1)
             if strupcase(cc) eq 'Y' then write = 1L
             print
          endif else write = 1L
          
          if write then begin
             openw, lun1, arcfile, /get_lun
             printf, lun1, '# Arc lamp list. '
             print_struct, ainfo, ['FILE','OBJECT','TYPE'], title, outstr, /strings
             printf, lun1, '# '+title
             for k = 0L, n_elements(outstr[0,*])-1L do printf, lun1, '  ', outstr[*,k], format='(A2,3A)'
             free_lun, lun1
          endif

          wmapfile = 'wmaplist_'+outroot+'.txt'
          if file_test(wmapfile,/regular) and (not keyword_set(overwrite)) then begin
             print, 'Overwrite existing file '+wmapfile+' [Y/N]? ', format='(A,$ )'
             cc = get_kbrd(1)
             if strupcase(cc) eq 'Y' then write = 1L
             print
          endif else write = 1L
          
          if write then begin
             openw, lun1, wmapfile, /get_lun
             printf, lun1, '# Wavelength map list. '
             printf, lun1, '# '
             for k = 0L, narcs-1L do printf, lun1, 'wmap_'+strmid(ainfo[k].file,0,strpos(ainfo[k].file,'.fits'))+'.idlsave'
             free_lun, lun1
          endif

       endelse 
          
    endif 
       
    popd
    
return
end

