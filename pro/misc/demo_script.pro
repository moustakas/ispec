pro demo_script, doplot=doplot
    
    datapath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='demo')
    pushd, datapath
    
; ----------------------------------------------------------------------
; example code for initializing input text files - do not use
; ----------------------------------------------------------------------
    
    imake_keywordfile, 'PA', comments='slit position angle', values='90.0', $
      root='a.', datapath=datapath, headfile='header_keywords.dat'
    ispec_makelists, 'a.', 'demo', datapath=datapath, /overwrite
    
; ----------------------------------------------------------------------

    paramfile = 'ibatch_demo.txt'
    procfile = 'objlist_demo.txt'
    crsplitfile = 'crsplits_demo.txt'
    crfile = 'crlist_demo.txt'
    tracefile = 'tracelist_demo.txt'
    stdfile = 'stdlist_demo.txt'
    calibfile = 'caliblist_demo.txt'

    headfile = 'header_keywords_demo.dat'
    iheader_keywords, headfile, datapath=datapath, /update, silent=silent

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PART 1 - using IBATCH
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; ----------------------------------------------------------------------
; [1] generate the flat field and the illumination correction
; ----------------------------------------------------------------------

    ibatch, paramfile, /makeflat, checkillum=checkillum, doplot=doplot

; ----------------------------------------------------------------------
; [2] run all the preliminary reductions on the data.  verify that
; OBJLIST_DEMO.TXT has no bias frames, dome flats, or darks 
; ----------------------------------------------------------------------

    readcol, procfile, proclist, format='A', /silent, comment='#'
    ibatch, paramfile, proclist=proclist, /ccdproc, checkoverscan=checkoverscan
    
; ----------------------------------------------------------------------
; [3] model the spatial distortion with the standard stars
; ----------------------------------------------------------------------

    readcol, tracefile, tracelist, format='A', /silent, comment='#'
    ibatch, paramfile, tracelist=tracelist, /distortion, doplot=doplot

; ----------------------------------------------------------------------
; [4] generate the wavelength solution
; ----------------------------------------------------------------------

    ibatch, paramfile, /arcfit, doplot=doplot

; ----------------------------------------------------------------------
; [5] combine cr-splits
; ----------------------------------------------------------------------

    ibatch, paramfile, crsplitfile=crsplitfile, /crsplits
    
; ----------------------------------------------------------------------
; [6] identify cosmic rays    
; ----------------------------------------------------------------------

    readcol, crfile, crlist, format='A', /silent, comment='#'
    ibatch, paramfile, crlist=crlist, /crclean

; ----------------------------------------------------------------------
; [7] generate the sensitivity function with a grey shift
; ----------------------------------------------------------------------

    readcol, stdfile, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile, stdlist=stdlist, /makesens, /grey, doplot=doplot

; ----------------------------------------------------------------------
; [8] calibrate all the data
; ----------------------------------------------------------------------

    readcol, calibfile, caliblist, format='A', /silent, comment='#'
    ibatch, paramfile, caliblist=caliblist, /calibrate

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PART 2 - using IALL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
; ----------------------------------------------------------------------
; do [1-8] in one call
; ----------------------------------------------------------------------
    
    iall, paramfile, procfile=procfile, tracefile=tracefile, crfile=crfile, $
      crsplitfile=crsplitfile, stdfile=stdfile, calibfile=calibfile, /ccdproc, $
      /crclean, /grey, /sens, /calibrate, doplot=doplot

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PART 3 - other handy routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
; ----------------------------------------------------------------------
; generate the data reduction web page
; ----------------------------------------------------------------------

    kgwebpage, 'demo', red_path='./', html_path=html_path, rootname='fdw', $
      /doskysub

; ----------------------------------------------------------------------
; compute the seeing
; ----------------------------------------------------------------------

    iseeing, 'seeing_demo.txt', title='Demo Seeing', psname='seeing_demo.ps', $
      /postscript, doplot=doplot
    
    popd
    
return
end
