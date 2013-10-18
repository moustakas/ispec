pro slit_function, flatname, datapath=datapath, trim=trim, postscript=postscript
;+
; NAME:
;	SLIT_FUNCTION
;
; PURPOSE:
;	Derive the slit function from a dome flat.
;
; CALLING SEQUENCE:
;	slit_function, flatname, [datapath=], postscript=postscript
;
; INPUTS:
;	flatname  - FITS file name of the dome flat
;
; OPTIONAL INPUTS:
;	datapath  - path to the dome flat (default PWD)
;	trim   - [x1,x2,y1,y2] defining the good data region
;                (zero-indexed)
;
; KEYWORD PARAMETERS:
;	postscript - generate postscript output
;
; PROCEDURES USED:
;	STRN(), SXPAR(), CWD()
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 October 22, U of A
;-

    if n_elements(flatname) eq 0L then begin
       print, 'Syntax - slit_function, flatname, [datapath=]'
       return
    endif
    
    if not keyword_set(datapath) then datapath = cwd()

    dome = readfits(flatname,dhead,/silent)
    dsize = size(dome,/dimension)

    if keyword_set(trim) then begin

       xy = trim

       xy[0] = xy[0] > 0L
       xy[1] = xy[1] < (dsize[0]-1L)
       xy[2] = xy[2] > 0L
       xy[3] = xy[3] < (dsize[1]-1L)

    endif else begin

       xy = lonarr(4)
       
       xy[0] = 0L
       xy[1] = dsize[0]-1L
       xy[2] = 0L
       xy[3] = dsize[1]-1L

    endelse

    ncols = xy[1]-xy[0]+1L
    nrows = xy[3]-xy[2]+1L    

    rowaxis = findgen(nrows)
    if strupcase(strn(sxpar(dhead,'OBJECT'))) ne 'DOME FLAT' then begin
       print, 'The slit function must be computed from a dome flat field.'
       return
    endif
    
    dflat = fltarr(ncols,nrows)
    for k = 0L, ncols-1L do dflat[k,*] = dome[xy[0]+k,xy[2]:xy[3]]/mean(dome[xy[0]+k,xy[2]:xy[3]])

    slit = total(dflat,1)/ncols

    root = strmid(flatname,0,strpos(flatname,'.fits'))
    qaplotname = 'qaplot_sfunction_'+root
    
    if keyword_set(postscript) then begin
       ps_open, datapath+qaplotname, /ps_fonts
       device, /inches, /times
    endif else window, 0, xs=550, ys=550

    plot, rowaxis, slit, xsty=3, ysty=3, xtit='Row Number', ytit='Slit Function', $
      thick=2.0, charsize=1.6, charthick=2.0, tit=strlowcase(root), ps=10

    if keyword_set(postscript) then ps_close

return
end
