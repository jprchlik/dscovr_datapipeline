; dsc_advanced_kp_fit_v2
;
; This is a fork from dsc_advanced_kp_fit.pro (11/14/2016). This
; version is intended to include
;
; (1) 1 minute and high res analysis together
; (2) integrated uncertainty/error analysis
;
; This retires the "mode" and "refine" keywords and settings in the previous version
;
; Uses
;    Moved dfc_loadflowangle_lookup_table.pro and dfc_geteffarea.pro to subroutines 
;    Also copied files restored in those programs to local directory and changed restore
;    path in both dfc_* files.
;   .compile /crater/projects/dscovr/simulator/dfc_v4/dfc_loadflowangle_lookup_table.pro (deprecated
;             2017/02/17 J. Prchlik)
;   .compile /crater/projects/dscovr/simulator/dfc_v4/dfc_geteffarea.pro (depcreated 
;             2017/02/17 J. Prchlik)
;   .compile zptfc_to_xyzgse.pro
;    @'/crater/utilities/idl/mike/idlstartup'
;
; SINGLE MAXWELLIAN:
; This routine presently uses a fit to the tallest tree above
; an empirical f = 0.025 noise, and extending 2 points to the low-
; energy side of the peak, with dynamic adjustment of range to capture
; hotter peaks
;
; DOUBLE MAXWELLIAN
; This "advanced" routine is intended to be a fork for more
; sophisticated approaches. It presently uses the single-maxwellian
; as a guess and expands point selection to +/- 7 channels from the
; peak for a two-maxwellian fit.
;
;
; USAGE
; dsc_advanced_kp_fit, year, doy, version, show=show, hold=hold, $
;                         save=save, rezero=rezero, $
;                         averaging_length = averaging_length, $
;                         flybacks=flybacks, $
;                         neg_offset=neg_offset, verbose=verbose, $
;                         clobber = clobber
;
;
; ABERATION CORRECTION NOW IN PROGRAM HOWEVER REPROCESSING OF ALL TIME NOT COMPLETE (PRCHLIK, J 2017/03/16); 
;
;NOTE ON CLOBBER KEYWORD
;If the IDL save file already exits at the version you want do not set clobber keyword.
;This way the program sees the idl exists and creates the cdf file from that file at the newest version (rapid)
;If you set clobber you will overwrite the previous idl save file and create the cdf file at the current version (slow).
;If you do not set the clobber keyword but the file is not yet processed then setting clobber is irrelavent because 
;a new idl and cdf file will be created at the most recent version (slow).
;
;
;

; ---------------------------------------------------------------------
;
; first fit function 
;
; ---------------------------------------------------------------------
pro one_gaussian, x, a, y, pder

n = n_elements(a)
nx = N_ELEMENTS(x)
if a[2] ne 0.0 then begin
   Z = (X-A[1])/A[2]            ;GET Z
   EZ = EXP(-Z^2/2.)            ;GAUSSIAN PART
endif else begin
   z = REPLICATE(FIX(100, TYPE=SIZE(x,/TYPE)), nx)
   ez = z*0
endelse

Y = (A[0]>0)*EZ ; note that Y is going to be forced positive

IF N_PARAMS(0) LE 3 THEN RETURN ;NEED PARTIAL?

PDER = FLTARR(nx, n) ;YES, MAKE ARRAY.
PDER[*,0] = EZ       ;COMPUTE PARTIALS
if a[2] ne 0. then PDER[*,1] = (A[0]>0) * EZ * Z/A[2]
PDER[*,2] = PDER[*,1] * Z 

end


; ---------------------------------------------------------------------
;
; second fit function 
;
; ---------------------------------------------------------------------
pro two_gaussians, x, a, y, pder

n = n_elements(a)
nx = N_ELEMENTS(x)
if a[2] ne 0.0 then begin
   Z1 = (X-A[1])/A[2]            ;GET Z
   EZ1 = EXP(-Z1^2/2.)            ;GAUSSIAN PART
endif else begin
   z1 = REPLICATE(FIX(100, TYPE=SIZE(x,/TYPE)), nx)
   ez1 = z1*0
endelse
Y1 = (A[0]>0)*EZ1 ; note that Y is going to be forced positive

if a[5] ne 0.0 then begin
   Z2 = (X-A[4])/A[5]            ;GET Z
   EZ2 = EXP(-Z2^2/2.)            ;GAUSSIAN PART
endif else begin
   z2 = REPLICATE(FIX(100, TYPE=SIZE(x,/TYPE)), nx)
   ez2 = z2*0
endelse
Y2 = (A[3]>0)*EZ2 ; note that Y is going to be forced positive

Y=Y1+Y2

IF N_PARAMS(0) LE 3 THEN RETURN ;NEED PARTIAL?

PDER = FLTARR(nx, n) ;YES, MAKE ARRAY.
PDER[*,0] = EZ1       ;COMPUTE PARTIALS
if a[2] ne 0. then PDER[*,1] = (A[0]>0) * EZ1 * Z1/A[2]
PDER[*,2] = PDER[*,1] * Z1 

PDER[*,3] = EZ2       ;COMPUTE PARTIALS
if a[5] ne 0. then PDER[*,4] = (A[0]>0) * EZ2 * Z2/A[5]
PDER[*,5] = PDER[*,4] * Z2 

end


; ---------------------------------------------------------------------
; subroutine to restore idl sav file iwth load flow angle table
;
; ---------------------------------------------------------------------

pro dfc_loadflowangle_lookup_table, x, y, EW, NS, tablefile=tablefile

if not keyword_set(tablefile) then $
   tablefile = '/crater/observatories/dscovr/code/modules/DSCOVR_FLOWANGLE_LOOKUP_HIRES_V1.idl'
restore, tablefile


end

; ---------------------------------------------------------------------
;
; subroutine to calculate flow angles from plate signals, using 
; the correction factors if desired.
;
; ---------------------------------------------------------------------

pro ABC2pt, sA,sB,sC, phi, theta, corrected=corrected

; this is the default correction matrix
p = [0., 1., 0., 0., $
     0., 0., 1., 0., $
     0., 0., 0., 1.]

if keyword_set(corrected) then begin
   restore, '/crater/observatories/dscovr/code/modules/dsc_plate_corrections.idl'
   p = dsc_plate_corrections
endif

; rescale the plate signals
sa1 = p[0] + p[1]*sa + p[2]*sb + p[3]*sc
sb1 = p[0+4] + p[1+4]*sa + p[2+4]*sb + p[3+4]*sc
sc1 = p[0+8] + p[1+8]*sa + p[2+8]*sb + p[3+8]*sc

x_fc = (sb1 - sa1)/(sa1+sb1+sc1)
y_fc = (sa1+sb1-sc1)/(sa1+sb1+sc1)
dfc_loadflowangle_lookup_table, xx, yy, EW, NS
flip_x = 2.*(x_FC gt 0) - 1.
x_FC = abs(x_fc)
xr = xx[*]
yr = yy[*]
xinterp = interpol(findgen(n_elements(xr)), xr, x_FC)
yinterp = interpol(findgen(n_elements(yr)), yr, y_FC)
phi = interpolate(EW, xinterp, yinterp, cubic=(-0.5)) * flip_x ; phi FC angle
theta = interpolate(NS, xinterp, yinterp, cubic=(-0.5)) ; theta FC angle

end

; ---------------------------------------------------------------------
;
; Apply fills to unphysical results
;
; ---------------------------------------------------------------------

pro rangechecks, dat

invalid = where(dat.vx.data gt (-200) or $
                dat.vx.data lt (-1200) or $
                abs(dat.vy.data) gt 250 or $
                abs(dat.vz.data) gt 250 or $
                dat.w.data lt 6. or $
                dat.w.data gt 160., ninval)

if ninval lt 1 then return
; (else)
for i = 2, n_tags(dat) - 2 do dat.(i).data[invalid] = dat.(i).fillval
for i = 2, n_tags(dat) - 2 do dat.(i).uncertainty[invalid] = dat.(i).fillval

end

;-------------------------------------------------
;
;USAGE
;running_med,x,y,medar,dmeda,npix=npix
;
;COMMENTS
;    Computes a running pixel median (default = 1, i.e. +/- 1 totaling 3 pixels)
;    
;-------------------------------------------------
pro running_med,y,medar,dmeda,npix=npix

if keyword_set(npix) then npix=npix  else npix=1   

sizer = n_elements(y)
medar = fltarr(sizer)

;good values to pull median from
find_val = where(y gt -9990.)
good_val = y[find_val]

for i=npix,sizer-npix-1 do begin
;get index in good_val to find median
    j = where(find_val eq i,gcnt)
;if value is already filled don't give user check flag
    if gcnt eq 0 then begin 
        medar[i] = y[i]
        continue
    endif
;get the median value near pixel
    ind_min = fix(j)-npix
    ind_max = fix(j)+npix
;prevent overreaching in array 
    while ind_max gt n_elements(good_val)-1 do ind_max = ind_max-1
;prevent under 0 in array 
    while ind_min lt 0 do ind_min = ind_min+1

    medvl = good_val[ind_min:ind_max]
;Count the number of acceptable pixels and where they are located to apply to case
    use = where(medvl gt -9990.,cnt)
    case 1 of
        (cnt gt 1): medar[i] = median(medvl[use]) ;find median of neighboring points
        (cnt eq 1): medar[i] = medvl[use];if only 1 good point assume it is the median
        (cnt eq 0): medar[i] = medvl[1] ;if they are all fill put the median as the fill val
    endcase
endfor



;correct for day ends
medar[0:npix-1] = medar[npix]
medar[sizer-npix:sizer-1] = medar[sizer-npix-1]

;find the difference between the median array and the measured value
dmeda = y-medar


end

;--------------------------------------------------
;
;USAGE
;replace = sig_replace,x,y,npix=npix,sigcut=sigcut,tol=tol
;
;COMMENTS
;Sends replacement array for values which signficantly differ from the running median
;--------------------------------------------------
function sig_replace,y,medar=medar,dmeda=dmeda,npix=npix,sigcut=sigcut,tol=tol
if keyword_set(npix) then npix = npix else npix = 1
if keyword_set(sigcut) then sigcut=sigcut else sigcut = 200.00;0.75
if keyword_set(tol) then tol = tol else tol=1 ; number of "bad" values allowed in a row


running_med,y,medar,dmeda,npix=npix



;Find values more than the sigcut number of sigma away
use = where((y gt -9990))
if n_elements(size(use)) lt 4 then return,-9999.0 ; end if no good spectrum data exists for day
;Tried difference sigma estimates; however none worked well so I am using a base 100km/s rejection
;sigmed = stddev(dmeda[use],/nan);/sqrt(n_elements(use))
sigmed = abs(dmeda[use]-median(dmeda[use]))
sigmed = median(sigmed)
replace = where((abs(dmeda) gt 100.) and (y gt -9990))

;plot,abs(dmeda)/sigmed,yrange=[0,550],ystyle=1,psym=6
;wait,20

;Check replace for any neighbor bad pixels if true assumed the "bad" assumption is false
;report only uniq values in rebad (replace bad array)
case 1 of
    ((n_elements(size(replace)) gt 3) and (n_elements(replace) ge 2)):begin
        check1 = replace[1:n_elements(replace)-1]-replace[0:n_elements(replace)-2]
        rebad1 = where(check1 gt 1,crebad1)
        if crebad1 gt 0 then begin 
            rebad1 = [rebad1,rebad1+1]
            srebad1= sort(rebad1)
            rebad1 = rebad1[srebad1]
            urebad1= uniq(rebad1)
            rebad  = rebad1[urebad1]
        endif
       
        if n_elements(replace) gt 2 then begin
            check2 = replace[2:n_elements(replace)-1]-replace[0:n_elements(replace)-3]
            rebad2 = where(check2 gt 2,crebad2)

            case 1 of 
                ((crebad1 gt 0) and (crebad2 gt 0)): begin
                    rebad2 = [rebad2,rebad2+2]
                    srebad2= sort(rebad2)
                    rebad2 = rebad2[srebad2]
                    urebad2= uniq(rebad2)
                    rebad2 = rebad2[urebad2]
                    rebad  = [rebad,rebad2]
                    srebad = sort(rebad)
                    rebad  = rebad[srebad]
                    urebad = uniq(rebad)
                    rebad  = rebad[urebad]
                end
                ((crebad1 eq 0) and (crebad2 gt 0)): begin
                    rebad2 = [rebad2,rebad2+2]
                    srebad2= sort(rebad2)
                    rebad2 = rebad2[srebad2]
                    urebad2= uniq(rebad2)
                    rebad2 = rebad2[urebad2]
                    rebad  = rebad2
                end
                ((crebad1 gt 0) and (crebad2 eq 0)): rebad = rebad
                ((crebad1 eq 0) and (crebad2 eq 0)): rebad = -9999
              
            endcase
        endif 
 

        ;return the lonely bad points
        if n_elements(rebad) gt 3 then replace = replace[rebad] else replace = -9999
    end 
    (n_elements(replace) eq 1): replace = replace
    else: replace = -9999
endcase
return,replace
end


; ---------------------------------------------------------------------
;
; SUBROUTINE that takes fit data output and 
; packages it into a nice data structure for 
; human use
;
; ---------------------------------------------------------------------
pro fit_output_structure, fitpars, fitsig, phi, theta, dphi, dtheta, year, spec_jd, result, notes=notes, $
                          fillval=fillval

; BASIC FITPAR EXTRACTION AND PREP: Full resolution
mp = 1.67262178d-27             ; SI units
kb = 1.3806488d-23              ;  SI units

u = reform(fitpars[1, *])
w = reform(sqrt(2.*fitpars[2, *]^2))
n = reform(fitpars[0,*])*w*sqrt(!pi)
T = 0.5*(mp/kb)*(1d3*w)^2
day = spec_jd - julday(1, 1, year, 0, 0, 0) + 1.
zptfc_to_xyzgse, [[u], [phi], [theta]], ugse, year, day

du = reform(fitsig[1, *])
dw = reform(fitsig[2, *])*sqrt(2.)
dn = sqrt(!pi)*sqrt((w*reform(fitsig[0, *]))^2 + (reform(fitpars[0, *])*dw)^2)
dT = 2.*T*dw/w

; the off-radial velocity components are high uncertainty
; we're going to treat them the same for the sake of uncertainty
; estimates
;fix for solar wind's aberration in Y component
solab = 29.78 ;km/s
ux = reform(ugse[*, 0])
uy = reform(ugse[*, 1])-solab ; Add (2017/04/07 Prchlik. J)
uz = reform(ugse[*, 2])
umag = sqrt(ux^2 + uy^2 + uz^2)
uperp = sqrt(umag^2 - ux^2)

alpha = sqrt(phi^2 + theta^2)
dalpha = sqrt(dphi^2 + dtheta^2)
dupar = du
duperp1 =  sqrt( (du^2)*(tan(phi*!dtor)^2) + $ 
               (ux^2)*((!dtor*dphi)^2)/(cos(phi*!dtor)^4))
duperp2 =  sqrt( (du^2)*(tan(theta*!dtor)^2) + $ 
               (ux^2)*((!dtor*dtheta)^2)/(cos(theta*!dtor)^4))
duperp = sqrt(duperp1^2 + duperp2^2)

dux = dupar
duy = duperp
duz = duperp

; fill bad data 
if not keyword_set(fillval) then fillval = -9999.
fill = where(fitpars[1, *] lt -3000.0, nfill)
if nfill gt 0 then begin
   u[fill] = fillval
   w[fill] = fillval
   n[fill] = fillval
   T[fill] = fillval
   du[fill] = fillval
   dw[fill] = fillval
   dn[fill] = fillval
   dT[fill] = fillval
   ux[fill] = fillval
   dux[fill] = fillval
   uy[fill] = fillval
   duy[fill] = fillval
   uz[fill] = fillval
   duz[fill] = fillval
endif

;fill data more than sigcut sigmas from the running median (default =6)
replace = sig_replace(ux)
if n_elements(size(replace)) gt 3 then begin
   u[replace] = fillval
   w[replace] = fillval
   n[replace] = fillval
   T[replace] = fillval
   du[replace] = fillval
   dw[replace] = fillval
   dn[replace] = fillval
   dT[replace] = fillval
   ux[replace] = fillval
   dux[replace] = fillval
   uy[replace] = fillval
   duy[replace] = fillval
   uz[replace] = fillval
   duz[replace] = fillval
endif


;fill = -9999.   ; come back to bad data flagging later or write a
;separate routine
;bad = where(ux gt 0 or n lt 0 or n gt 100 or T lt 1e3 or T gt 1e7)
if not keyword_set(notes) then notes =  ['Prepared with dsc_advanced_kp_fit_v2', $
                   'mstevens@cfa.harvard.edu '+systime()] 
result = {year:year, time:{data:day, name: 'time stamp', $
                           units:'fractional day of year [UT]'}, $
           vx:{data:ux, name: 'proton velocity x-component', $
               units:'km/s', frame:'GSE', uncertainty:dux, fillval:fillval}, $
           vy:{data:uy, name: 'proton velocity y-component', $
               units:'km/s', frame:'GSE', uncertainty:duy, fillval:fillval}, $
           vz:{data:uz, name: 'proton velocity z-component', $
               units:'km/s', frame:'GSE', uncertainty:duz, fillval:fillval}, $
           n: {data:n, name: 'proton density', units:'cm^-3', $
              uncertainty:dn, fillval:fillval}, $
           w: {data:w, name: 'proton most probable thermal speed', units:'km/s', $
              uncertainty:dw, fillval:fillval}, $
           T: {data:T, name: 'proton temperature', units:'Kelvin', $
              uncertainty:dT, fillval:fillval}, $
           notes: notes     }

end


; ---------------------------------------------------------------------
;
; SUBROUTINE for decimating  oversampled smoothed data onto 1 minute grid
;
; ---------------------------------------------------------------------
pro regrid_1minute, invars, outvars, theday, doy, outdoy=outdoy, fillval=fillval

; We need to filter out bad points, and resample onto an appropriate grid
 num_minutes = 24.*60.
 minutes_grid = lindgen(num_minutes) ; this is minutes from 00:00 on the first day of the selection
 doy_grid = min(fix(theday)) + minutes_grid/(24.*60.)
 num_vars = n_elements(invars[0, *])

 outvars = fltarr(num_minutes, num_vars)

 for i = 0, num_vars - 1 do $
    outvars[*, i] = interpol(median(reform(invars[*, i]), 15), doy, doy_grid)
    
;
minutes = round((doy - min(round(doy)))*24.*60.)
unique_minutes = minutes[uniq(minutes)]

; Here's a little differencing algorithm
a = minutes_grid
b = unique_minutes
mina = Min(a, Max=maxa)
minb = Min(b, Max=maxb)
IF (minb GT maxa) OR (maxb LT mina) THEN diff= a ;No intersection...
r = Where((Histogram(a, Min=mina, Max=maxa) NE 0) AND $
          (Histogram(b, Min=mina, Max=maxa) EQ 0), count)
IF count eq 0 THEN diff =-1 ELSE diff= r + mina

; r contains the values in minutes_grid that aren't present in the 
; dscovr unique minutes. These are the ones we need to black out.
grid_tofill = r
if not keyword_set(fillval) then fillval = -9999.99
if r[0] ne (-1) then outvars[r, *] = fillval
outdoy = doy_grid

end

; ---------------------------------------------------------------------
;
; SUBROUTINE for loading current spectrum data
;
; ---------------------------------------------------------------------
pro get_spectra, year, doy, ia, ib, ic, spec_jd, verbose=verbose

; set up desired date range
jd = doy + julday(1, 1, year, 0, 0, 0) - 1
jdb4 = jd - 1
caldat, jdb4, mb4, db4, yb4
doyb4 = jdb4 - julday(1, 1, yb4, 0, 0, 0) + 1.
jdmin = jd - 0.007              ; 10 minute margin for interpolations
jdmax = jd + 1.007              ; 10 minute margin for interpolations

; (1) INITIALIZATIONS AND LOADS
; We may need to also load a bit from the day before
restore, '/crater/observatories/dscovr/code/modules/dsc_Vwindows.idl'
loadct, 39
device, decomposed = 0
ddd = string(doy, format = '(I03)')
dddb4 = string(doyb4, format = '(I03)')
yyyy = string(year, format = '(I4)')
yyyyb4 = string(yb4, format = '(I4)')

path = '/crater/observatories/dscovr/plasmag/l2/idl/ionspec/' + yyyy + '/'
path2 =  '/crater/observatories/dscovr/plasmag/l2/idl/ionspec/' + yyyyb4 + '/'
specfile = 'dsc_fc_ionspec_'+yyyy+'_'+ddd+'.idl'
specfile2 =  'dsc_fc_ionspec_'+yyyyb4+'_'+dddb4+'.idl'
result = file_search(path+specfile, count = count)
resultb4 = file_search(path2+specfile2, count = countb4)

if keyword_set(verbose) then begin
   print, '     dsc_advanced_kp_fit:get_spectra:  ' 
   print, '     Data files found, loading:'
   print, '          '+resultb4
   print, '          '+result
endif


; Handle cases where the target file is not found
if count eq 0 then begin
   print, '------------------------------------------------------------------'
   print, '     dsc_advanced_kp_fit:get_spectra:  ' 
   print, '          Error, no ion spec file found for ' + ddd + '/' + yyyy
   print, '          Checking for data in previous day file'
   if countb4 eq 0 then begin
      print, '          No data in previous day, either. returning'
      print, '------------------------------------------------------------------'
   endif
endif else begin
   restore, result
   tk = where(spec_jd gt jdmin and spec_jd lt jdmax, ntk)
   if ntk gt 0 then begin
      ia = spec[tk].ia
      ib = spec[tk].ib
      ic = spec[tk].ic
      hold_specjd = spec_jd[tk]
   endif
endelse

; Handle data from the previous day file that 
; dribbles into the target day
if countb4 gt 0 then begin
   restore, resultb4
   tk = where(spec_jd gt jdmin and spec_jd lt jdmax, ntk)
   if ntk gt 0 then begin
      iab4 = spec[tk].ia
      ibb4 = spec[tk].ib
      icb4 = spec[tk].ic
      hold_specjdb4 = spec_jd[tk]
      if count gt 0 then begin
         ia=[[iab4], [ia]]
         ib=[[ibb4], [ib]]
         ic=[[icb4], [ic]]
         hold_specjd = [hold_specjdb4, hold_specjd]
      endif else begin
         ia=iab4
         ib=ibb4
         ic=icb4
         hold_specjd = hold_specjdb4
      endelse
   endif
endif

spec_jd = hold_specjd
 
end

; ---------------------------------------------------------------------
;
; SUBROUTINE for saving the outputs
;
; ---------------------------------------------------------------------
pro advkp_save, adv_kp, dfc_kp, dfc_kp_1min, yyyy, ddd

; set up output file names
   outpath =  '/crater/observatories/dscovr/plasmag/l2/idl/advkp/' + yyyy + '/' 
   prefix = 'dsc_fc_advkp'
   outfile = prefix + '_' + yyyy + '_' + ddd+'.idl'

   outpath_highres =  '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/fullres/' 
   prefix_highres = 'dsc_fc_advkp_fullres'
   outfile_highres = prefix_highres + '_' + yyyy + '_' + ddd+'.idl'

   outpath_1min =  '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute/' 
   prefix_1min = 'dsc_fc_advkp_1minute'
   outfile_1min = prefix_1min + '_' + yyyy + '_' + ddd+'.idl'

   save, dfc_kp, filename = outpath_highres+outfile_highres
   save, dfc_kp_1min, filename = outpath_1min+outfile_1min
   save, adv_kp, filename = outpath+outfile
end


; ---------------------------------------------------------------------
;
; SUBROUTINE for obtaining the effective area of the FC
;
; ---------------------------------------------------------------------
; DSCOVR FARADAY CUP SIMULATION ALGORITHMS
; Michael L. Stevens
; mstevens@cfa.harvard.edu

; load an array of dscovr FC effective area, in cc/km, as a function of 
; incidence angle. The returned array is indexed in tenths of a degree, i.e.
; 
; units are cm^3 per km

function dfc_geteffarea, phi_deg, theta_deg, area_file=area_file, show=show, thetas=thetas, phis=phis, effarea_dscovrsim=effarea_dscovrsim

if (min(phi_deg^2 + theta_deg^2) ge 65.^2) then return, 0.*phi_deg

if not keyword_set(area_file) then area_file = '/crater/observatories/dscovr/code/modules/DSCOVR_FC_EFFAREA_LOOKUP_V1.idl'

if not keyword_set(thetas) or (not keyword_set(phis)) or (not keyword_set(effarea_dscovrsim)) then restore, area_file
; table is symmetrical in phi
absphi_deg = abs(phi_deg)

; interpolate on the table for the flow angle at the lookup values
phir = phis[*, 0]
thetar = thetas[0, *]
phiinterp = interpol(findgen(n_elements(phir)), phir, absphi_deg)
thetainterp = interpol(findgen(n_elements(thetar)), thetar, theta_deg)

effarea = interpolate(effarea_dscovrsim, phiinterp, thetainterp, /cubic)

if keyword_set(show) then begin
  levels = [0.99*effarea, 1.01*effarea]
  contour, effarea_dscovrsim, phis, thetas, levels = levels, /follow, c_linestyle = 1
  contour, effarea_dscovrsim, phis, thetas, levels = [effarea], /follow, /overplot
  oplot, [absphi_deg], [theta_deg], psym = 2, thick = 2, symsize = 2
endif  

return, effarea

end



pro dfc_geteffarea_windsim, eff_area

r_cm = 4.0; centimeters
a_cm2 = !dpi*r_cm^2
grid_transp = 0.9025
ngrids = 8.
A0_cm2 = a_cm2 * (grid_transp^ngrids)
A0_cm3km1 = A0_cm2*1d5

alpha = 0.1*dindgen(900)
eff_area = dblarr(900)

; As of 10-MAR-2014, we have not identified any differences between the dscovr cup
; and the Wind cups that would affect the effective area except for the geometry of
; the collector. The working lookup table for the Wind cups implies a [2D] grid 
; transparency of 95.7%, which is higher than we thought [expected 90%], but we have 
; no measurements at present to imply otherwise.
;
; This needs to be corrected because dscovr has 8 grids to Wind's 9
;
;
; Because the dscovr grid transparencies are not known directly at present and because the 
; best known aperture radius is equal to that of Wind, we use the Wind effective 
; area table directly for dscovr until such time as we find reason to do otherwise.
;
eff_area[  0]=3.3820000e+06&eff_area[  1]=3.3821000e+06&eff_area[  2]=3.3822000e+06
eff_area[  3]=3.3823000e+06&eff_area[  4]=3.3824000e+06&eff_area[  5]=3.3825000e+06
eff_area[  6]=3.3826000e+06&eff_area[  7]=3.3827000e+06&eff_area[  8]=3.3828000e+06
eff_area[  9]=3.3829000e+06&eff_area[ 10]=3.3830000e+06&eff_area[ 11]=3.3830000e+06
eff_area[ 12]=3.3830000e+06&eff_area[ 13]=3.3830000e+06&eff_area[ 14]=3.3830000e+06
eff_area[ 15]=3.3830000e+06&eff_area[ 16]=3.3830000e+06&eff_area[ 17]=3.3830000e+06
eff_area[ 18]=3.3830000e+06&eff_area[ 19]=3.3830000e+06&eff_area[ 20]=3.3830000e+06
eff_area[ 21]=3.3829000e+06&eff_area[ 22]=3.3828000e+06&eff_area[ 23]=3.3827000e+06
eff_area[ 24]=3.3826000e+06&eff_area[ 25]=3.3825000e+06&eff_area[ 26]=3.3824000e+06
eff_area[ 27]=3.3823000e+06&eff_area[ 28]=3.3822000e+06&eff_area[ 29]=3.3821000e+06
eff_area[ 30]=3.3820000e+06&eff_area[ 31]=3.3819000e+06&eff_area[ 32]=3.3818000e+06
eff_area[ 33]=3.3817000e+06&eff_area[ 34]=3.3816000e+06&eff_area[ 35]=3.3815000e+06
eff_area[ 36]=3.3814000e+06&eff_area[ 37]=3.3813000e+06&eff_area[ 38]=3.3812000e+06
eff_area[ 39]=3.3811000e+06&eff_area[ 40]=3.3810000e+06&eff_area[ 41]=3.3809000e+06
eff_area[ 42]=3.3808000e+06&eff_area[ 43]=3.3807000e+06&eff_area[ 44]=3.3806000e+06
eff_area[ 45]=3.3805000e+06&eff_area[ 46]=3.3804000e+06&eff_area[ 47]=3.3803000e+06
eff_area[ 48]=3.3802000e+06&eff_area[ 49]=3.3801000e+06&eff_area[ 50]=3.3800000e+06
eff_area[ 51]=3.3798000e+06&eff_area[ 52]=3.3796000e+06&eff_area[ 53]=3.3794000e+06
eff_area[ 54]=3.3792000e+06&eff_area[ 55]=3.3790000e+06&eff_area[ 56]=3.3788000e+06
eff_area[ 57]=3.3786000e+06&eff_area[ 58]=3.3784000e+06&eff_area[ 59]=3.3782000e+06
eff_area[ 60]=3.3780000e+06&eff_area[ 61]=3.3779000e+06&eff_area[ 62]=3.3778000e+06
eff_area[ 63]=3.3777000e+06&eff_area[ 64]=3.3776000e+06&eff_area[ 65]=3.3775000e+06
eff_area[ 66]=3.3774000e+06&eff_area[ 67]=3.3773000e+06&eff_area[ 68]=3.3772000e+06
eff_area[ 69]=3.3771000e+06&eff_area[ 70]=3.3770000e+06&eff_area[ 71]=3.3769000e+06
eff_area[ 72]=3.3768000e+06&eff_area[ 73]=3.3767000e+06&eff_area[ 74]=3.3766000e+06
eff_area[ 75]=3.3765000e+06&eff_area[ 76]=3.3764000e+06&eff_area[ 77]=3.3763000e+06
eff_area[ 78]=3.3762000e+06&eff_area[ 79]=3.3761000e+06&eff_area[ 80]=3.3760000e+06
eff_area[ 81]=3.3758000e+06&eff_area[ 82]=3.3756000e+06&eff_area[ 83]=3.3754000e+06
eff_area[ 84]=3.3752000e+06&eff_area[ 85]=3.3750000e+06&eff_area[ 86]=3.3748000e+06
eff_area[ 87]=3.3746000e+06&eff_area[ 88]=3.3744000e+06&eff_area[ 89]=3.3742000e+06
eff_area[ 90]=3.3740000e+06&eff_area[ 91]=3.3738000e+06&eff_area[ 92]=3.3736000e+06
eff_area[ 93]=3.3734000e+06&eff_area[ 94]=3.3732000e+06&eff_area[ 95]=3.3730000e+06
eff_area[ 96]=3.3728000e+06&eff_area[ 97]=3.3726000e+06&eff_area[ 98]=3.3724000e+06
eff_area[ 99]=3.3722000e+06&eff_area[100]=3.3720000e+06&eff_area[101]=3.3717000e+06
eff_area[102]=3.3714000e+06&eff_area[103]=3.3711000e+06&eff_area[104]=3.3708000e+06
eff_area[105]=3.3705000e+06&eff_area[106]=3.3702000e+06&eff_area[107]=3.3699000e+06
eff_area[108]=3.3696000e+06&eff_area[109]=3.3693000e+06&eff_area[110]=3.3690000e+06
eff_area[111]=3.3689000e+06&eff_area[112]=3.3688000e+06&eff_area[113]=3.3687000e+06
eff_area[114]=3.3686000e+06&eff_area[115]=3.3685000e+06&eff_area[116]=3.3684000e+06
eff_area[117]=3.3683000e+06&eff_area[118]=3.3682000e+06&eff_area[119]=3.3681000e+06
eff_area[120]=3.3680000e+06&eff_area[121]=3.3676000e+06&eff_area[122]=3.3672000e+06
eff_area[123]=3.3668000e+06&eff_area[124]=3.3664000e+06&eff_area[125]=3.3660000e+06
eff_area[126]=3.3656000e+06&eff_area[127]=3.3652000e+06&eff_area[128]=3.3648000e+06
eff_area[129]=3.3644000e+06&eff_area[130]=3.3640000e+06&eff_area[131]=3.3638000e+06
eff_area[132]=3.3636000e+06&eff_area[133]=3.3634000e+06&eff_area[134]=3.3632000e+06
eff_area[135]=3.3630000e+06&eff_area[136]=3.3628000e+06&eff_area[137]=3.3626000e+06
eff_area[138]=3.3624000e+06&eff_area[139]=3.3622000e+06&eff_area[140]=3.3620000e+06
eff_area[141]=3.3617000e+06&eff_area[142]=3.3614000e+06&eff_area[143]=3.3611000e+06
eff_area[144]=3.3608000e+06&eff_area[145]=3.3605000e+06&eff_area[146]=3.3602000e+06
eff_area[147]=3.3599000e+06&eff_area[148]=3.3596000e+06&eff_area[149]=3.3593000e+06
eff_area[150]=3.3590000e+06&eff_area[151]=3.3586000e+06&eff_area[152]=3.3582000e+06
eff_area[153]=3.3578000e+06&eff_area[154]=3.3574000e+06&eff_area[155]=3.3570000e+06
eff_area[156]=3.3566000e+06&eff_area[157]=3.3562000e+06&eff_area[158]=3.3558000e+06
eff_area[159]=3.3554000e+06&eff_area[160]=3.3550000e+06&eff_area[161]=3.3546000e+06
eff_area[162]=3.3542000e+06&eff_area[163]=3.3538000e+06&eff_area[164]=3.3534000e+06
eff_area[165]=3.3530000e+06&eff_area[166]=3.3526000e+06&eff_area[167]=3.3522000e+06
eff_area[168]=3.3518000e+06&eff_area[169]=3.3514000e+06&eff_area[170]=3.3510000e+06
eff_area[171]=3.3506000e+06&eff_area[172]=3.3502000e+06&eff_area[173]=3.3498000e+06
eff_area[174]=3.3494000e+06&eff_area[175]=3.3490000e+06&eff_area[176]=3.3486000e+06
eff_area[177]=3.3482000e+06&eff_area[178]=3.3478000e+06&eff_area[179]=3.3474000e+06
eff_area[180]=3.3470000e+06&eff_area[181]=3.3466000e+06&eff_area[182]=3.3462000e+06
eff_area[183]=3.3458000e+06&eff_area[184]=3.3454000e+06&eff_area[185]=3.3450000e+06
eff_area[186]=3.3446000e+06&eff_area[187]=3.3442000e+06&eff_area[188]=3.3438000e+06
eff_area[189]=3.3434000e+06&eff_area[190]=3.3430000e+06&eff_area[191]=3.3425700e+06
eff_area[192]=3.3421400e+06&eff_area[193]=3.3417100e+06&eff_area[194]=3.3412800e+06
eff_area[195]=3.3408500e+06&eff_area[196]=3.3404200e+06&eff_area[197]=3.3399900e+06
eff_area[198]=3.3395600e+06&eff_area[199]=3.3391300e+06&eff_area[200]=3.3387000e+06
eff_area[201]=3.3382400e+06&eff_area[202]=3.3377800e+06&eff_area[203]=3.3373200e+06
eff_area[204]=3.3368600e+06&eff_area[205]=3.3364000e+06&eff_area[206]=3.3359400e+06
eff_area[207]=3.3354800e+06&eff_area[208]=3.3350200e+06&eff_area[209]=3.3345600e+06
eff_area[210]=3.3341000e+06&eff_area[211]=3.3336200e+06&eff_area[212]=3.3331400e+06
eff_area[213]=3.3326600e+06&eff_area[214]=3.3321800e+06&eff_area[215]=3.3317000e+06
eff_area[216]=3.3312200e+06&eff_area[217]=3.3307400e+06&eff_area[218]=3.3302600e+06
eff_area[219]=3.3297800e+06&eff_area[220]=3.3293000e+06&eff_area[221]=3.3288000e+06
eff_area[222]=3.3283000e+06&eff_area[223]=3.3278000e+06&eff_area[224]=3.3273000e+06
eff_area[225]=3.3268000e+06&eff_area[226]=3.3263000e+06&eff_area[227]=3.3258000e+06
eff_area[228]=3.3253000e+06&eff_area[229]=3.3248000e+06&eff_area[230]=3.3243000e+06
eff_area[231]=3.3236900e+06&eff_area[232]=3.3230800e+06&eff_area[233]=3.3224700e+06
eff_area[234]=3.3218600e+06&eff_area[235]=3.3212500e+06&eff_area[236]=3.3206400e+06
eff_area[237]=3.3200300e+06&eff_area[238]=3.3194200e+06&eff_area[239]=3.3188100e+06
eff_area[240]=3.3182000e+06&eff_area[241]=3.3176600e+06&eff_area[242]=3.3171200e+06
eff_area[243]=3.3165800e+06&eff_area[244]=3.3160400e+06&eff_area[245]=3.3155000e+06
eff_area[246]=3.3149600e+06&eff_area[247]=3.3144200e+06&eff_area[248]=3.3138800e+06
eff_area[249]=3.3133400e+06&eff_area[250]=3.3128000e+06&eff_area[251]=3.3121500e+06
eff_area[252]=3.3115000e+06&eff_area[253]=3.3108500e+06&eff_area[254]=3.3102000e+06
eff_area[255]=3.3095500e+06&eff_area[256]=3.3089000e+06&eff_area[257]=3.3082500e+06
eff_area[258]=3.3076000e+06&eff_area[259]=3.3069500e+06&eff_area[260]=3.3063000e+06
eff_area[261]=3.3056300e+06&eff_area[262]=3.3049600e+06&eff_area[263]=3.3042900e+06
eff_area[264]=3.3036200e+06&eff_area[265]=3.3029500e+06&eff_area[266]=3.3022800e+06
eff_area[267]=3.3016100e+06&eff_area[268]=3.3009400e+06&eff_area[269]=3.3002700e+06
eff_area[270]=3.2996000e+06&eff_area[271]=3.2989200e+06&eff_area[272]=3.2982400e+06
eff_area[273]=3.2975600e+06&eff_area[274]=3.2968800e+06&eff_area[275]=3.2962000e+06
eff_area[276]=3.2955200e+06&eff_area[277]=3.2948400e+06&eff_area[278]=3.2941600e+06
eff_area[279]=3.2934800e+06&eff_area[280]=3.2928000e+06&eff_area[281]=3.2921100e+06
eff_area[282]=3.2914200e+06&eff_area[283]=3.2907300e+06&eff_area[284]=3.2900400e+06
eff_area[285]=3.2893500e+06&eff_area[286]=3.2886600e+06&eff_area[287]=3.2879700e+06
eff_area[288]=3.2872800e+06&eff_area[289]=3.2865900e+06&eff_area[290]=3.2859000e+06
eff_area[291]=3.2850900e+06&eff_area[292]=3.2842800e+06&eff_area[293]=3.2834700e+06
eff_area[294]=3.2826600e+06&eff_area[295]=3.2818500e+06&eff_area[296]=3.2810400e+06
eff_area[297]=3.2802300e+06&eff_area[298]=3.2794200e+06&eff_area[299]=3.2786100e+06
eff_area[300]=3.2778000e+06&eff_area[301]=3.2770900e+06&eff_area[302]=3.2763800e+06
eff_area[303]=3.2756700e+06&eff_area[304]=3.2749600e+06&eff_area[305]=3.2742500e+06
eff_area[306]=3.2735400e+06&eff_area[307]=3.2728300e+06&eff_area[308]=3.2721200e+06
eff_area[309]=3.2714100e+06&eff_area[310]=3.2707000e+06&eff_area[311]=3.2697900e+06
eff_area[312]=3.2688800e+06&eff_area[313]=3.2679700e+06&eff_area[314]=3.2670600e+06
eff_area[315]=3.2661500e+06&eff_area[316]=3.2652400e+06&eff_area[317]=3.2643300e+06
eff_area[318]=3.2634200e+06&eff_area[319]=3.2625100e+06&eff_area[320]=3.2616000e+06
eff_area[321]=3.2607900e+06&eff_area[322]=3.2599800e+06&eff_area[323]=3.2591700e+06
eff_area[324]=3.2583600e+06&eff_area[325]=3.2575500e+06&eff_area[326]=3.2567400e+06
eff_area[327]=3.2559300e+06&eff_area[328]=3.2551200e+06&eff_area[329]=3.2543100e+06
eff_area[330]=3.2535000e+06&eff_area[331]=3.2526000e+06&eff_area[332]=3.2517000e+06
eff_area[333]=3.2508000e+06&eff_area[334]=3.2499000e+06&eff_area[335]=3.2490000e+06
eff_area[336]=3.2481000e+06&eff_area[337]=3.2472000e+06&eff_area[338]=3.2463000e+06
eff_area[339]=3.2454000e+06&eff_area[340]=3.2445000e+06&eff_area[341]=3.2435100e+06
eff_area[342]=3.2425200e+06&eff_area[343]=3.2415300e+06&eff_area[344]=3.2405400e+06
eff_area[345]=3.2395500e+06&eff_area[346]=3.2385600e+06&eff_area[347]=3.2375700e+06
eff_area[348]=3.2365800e+06&eff_area[349]=3.2355900e+06&eff_area[350]=3.2346000e+06
eff_area[351]=3.2336300e+06&eff_area[352]=3.2326600e+06&eff_area[353]=3.2316900e+06
eff_area[354]=3.2307200e+06&eff_area[355]=3.2297500e+06&eff_area[356]=3.2287800e+06
eff_area[357]=3.2278100e+06&eff_area[358]=3.2268400e+06&eff_area[359]=3.2258700e+06
eff_area[360]=3.2249000e+06&eff_area[361]=3.2224200e+06&eff_area[362]=3.2199400e+06
eff_area[363]=3.2174600e+06&eff_area[364]=3.2149800e+06&eff_area[365]=3.2125000e+06
eff_area[366]=3.2100200e+06&eff_area[367]=3.2075400e+06&eff_area[368]=3.2050600e+06
eff_area[369]=3.2025800e+06&eff_area[370]=3.2001000e+06&eff_area[371]=3.1962400e+06
eff_area[372]=3.1923800e+06&eff_area[373]=3.1885200e+06&eff_area[374]=3.1846600e+06
eff_area[375]=3.1808000e+06&eff_area[376]=3.1769400e+06&eff_area[377]=3.1730800e+06
eff_area[378]=3.1692200e+06&eff_area[379]=3.1653600e+06&eff_area[380]=3.1615000e+06
eff_area[381]=3.1567500e+06&eff_area[382]=3.1520000e+06&eff_area[383]=3.1472500e+06
eff_area[384]=3.1425000e+06&eff_area[385]=3.1377500e+06&eff_area[386]=3.1330000e+06
eff_area[387]=3.1282500e+06&eff_area[388]=3.1235000e+06&eff_area[389]=3.1187500e+06
eff_area[390]=3.1140000e+06&eff_area[391]=3.1084820e+06&eff_area[392]=3.1029640e+06
eff_area[393]=3.0974460e+06&eff_area[394]=3.0919280e+06&eff_area[395]=3.0864100e+06
eff_area[396]=3.0808920e+06&eff_area[397]=3.0753740e+06&eff_area[398]=3.0698560e+06
eff_area[399]=3.0643380e+06&eff_area[400]=3.0588200e+06&eff_area[401]=3.0526550e+06
eff_area[402]=3.0464900e+06&eff_area[403]=3.0403250e+06&eff_area[404]=3.0341600e+06
eff_area[405]=3.0279950e+06&eff_area[406]=3.0218300e+06&eff_area[407]=3.0156650e+06
eff_area[408]=3.0095000e+06&eff_area[409]=3.0033350e+06&eff_area[410]=2.9971700e+06
eff_area[411]=2.9904530e+06&eff_area[412]=2.9837360e+06&eff_area[413]=2.9770190e+06
eff_area[414]=2.9703020e+06&eff_area[415]=2.9635850e+06&eff_area[416]=2.9568680e+06
eff_area[417]=2.9501510e+06&eff_area[418]=2.9434340e+06&eff_area[419]=2.9367170e+06
eff_area[420]=2.9300000e+06&eff_area[421]=2.9227000e+06&eff_area[422]=2.9154000e+06
eff_area[423]=2.9081000e+06&eff_area[424]=2.9008000e+06&eff_area[425]=2.8935000e+06
eff_area[426]=2.8862000e+06&eff_area[427]=2.8789000e+06&eff_area[428]=2.8716000e+06
eff_area[429]=2.8643000e+06&eff_area[430]=2.8570000e+06&eff_area[431]=2.8492000e+06
eff_area[432]=2.8414000e+06&eff_area[433]=2.8336000e+06&eff_area[434]=2.8258000e+06
eff_area[435]=2.8180000e+06&eff_area[436]=2.8102000e+06&eff_area[437]=2.8024000e+06
eff_area[438]=2.7946000e+06&eff_area[439]=2.7868000e+06&eff_area[440]=2.7790000e+06
eff_area[441]=2.7705000e+06&eff_area[442]=2.7620000e+06&eff_area[443]=2.7535000e+06
eff_area[444]=2.7450000e+06&eff_area[445]=2.7365000e+06&eff_area[446]=2.7280000e+06
eff_area[447]=2.7195000e+06&eff_area[448]=2.7110000e+06&eff_area[449]=2.7025000e+06
eff_area[450]=2.6940000e+06&eff_area[451]=2.6833000e+06&eff_area[452]=2.6725999e+06
eff_area[453]=2.6618999e+06&eff_area[454]=2.6511999e+06&eff_area[455]=2.6404999e+06
eff_area[456]=2.6297998e+06&eff_area[457]=2.6190998e+06&eff_area[458]=2.6083998e+06
eff_area[459]=2.5976998e+06&eff_area[460]=2.5869998e+06&eff_area[461]=2.5747998e+06
eff_area[462]=2.5625998e+06&eff_area[463]=2.5503998e+06&eff_area[464]=2.5381999e+06
eff_area[465]=2.5259999e+06&eff_area[466]=2.5137999e+06&eff_area[467]=2.5015999e+06
eff_area[468]=2.4894000e+06&eff_area[469]=2.4772000e+06&eff_area[470]=2.4650000e+06
eff_area[471]=2.4514999e+06&eff_area[472]=2.4379999e+06&eff_area[473]=2.4244999e+06
eff_area[474]=2.4109998e+06&eff_area[475]=2.3974998e+06&eff_area[476]=2.3839997e+06
eff_area[477]=2.3704996e+06&eff_area[478]=2.3569996e+06&eff_area[479]=2.3434996e+06
eff_area[480]=2.3299995e+06&eff_area[481]=2.3152995e+06&eff_area[482]=2.3005996e+06
eff_area[483]=2.2858997e+06&eff_area[484]=2.2711997e+06&eff_area[485]=2.2564998e+06
eff_area[486]=2.2417998e+06&eff_area[487]=2.2270998e+06&eff_area[488]=2.2123999e+06
eff_area[489]=2.1977000e+06&eff_area[490]=2.1830000e+06&eff_area[491]=2.1673000e+06
eff_area[492]=2.1515999e+06&eff_area[493]=2.1358999e+06&eff_area[494]=2.1201999e+06
eff_area[495]=2.1044998e+06&eff_area[496]=2.0887998e+06&eff_area[497]=2.0730997e+06
eff_area[498]=2.0573997e+06&eff_area[499]=2.0416997e+06&eff_area[500]=2.0259996e+06
eff_area[501]=2.0092997e+06&eff_area[502]=1.9925997e+06&eff_area[503]=1.9758998e+06
eff_area[504]=1.9591998e+06&eff_area[505]=1.9424999e+06&eff_area[506]=1.9257999e+06
eff_area[507]=1.9091000e+06&eff_area[508]=1.8924000e+06&eff_area[509]=1.8757001e+06
eff_area[510]=1.8590001e+06&eff_area[511]=1.8414001e+06&eff_area[512]=1.8238000e+06
eff_area[513]=1.8062000e+06&eff_area[514]=1.7885999e+06&eff_area[515]=1.7709999e+06
eff_area[516]=1.7533998e+06&eff_area[517]=1.7357998e+06&eff_area[518]=1.7181997e+06
eff_area[519]=1.7005997e+06&eff_area[520]=1.6829996e+06&eff_area[521]=1.6643997e+06
eff_area[522]=1.6457997e+06&eff_area[523]=1.6271998e+06&eff_area[524]=1.6085998e+06
eff_area[525]=1.5899999e+06&eff_area[526]=1.5713999e+06&eff_area[527]=1.5528000e+06
eff_area[528]=1.5342000e+06&eff_area[529]=1.5156001e+06&eff_area[530]=1.4970001e+06
eff_area[531]=1.4775001e+06&eff_area[532]=1.4580000e+06&eff_area[533]=1.4385000e+06
eff_area[534]=1.4189999e+06&eff_area[535]=1.3994999e+06&eff_area[536]=1.3799998e+06
eff_area[537]=1.3604998e+06&eff_area[538]=1.3409997e+06&eff_area[539]=1.3214997e+06
eff_area[540]=1.3019996e+06&eff_area[541]=1.2816997e+06&eff_area[542]=1.2613997e+06
eff_area[543]=1.2410998e+06&eff_area[544]=1.2207998e+06&eff_area[545]=1.2004999e+06
eff_area[546]=1.1801999e+06&eff_area[547]=1.1599000e+06&eff_area[548]=1.1396000e+06
eff_area[549]=1.1193001e+06&eff_area[550]=1.0990001e+06&eff_area[551]=1.0778801e+06
eff_area[552]=1.0567600e+06&eff_area[553]=1.0356400e+06&eff_area[554]=1.0145199e+06
eff_area[555]=9.9339984e+05&eff_area[556]=9.7227979e+05&eff_area[557]=9.5115973e+05
eff_area[558]=9.3003968e+05&eff_area[559]=9.0891962e+05&eff_area[560]=8.8779956e+05
eff_area[561]=8.6586962e+05&eff_area[562]=8.4393969e+05&eff_area[563]=8.2200975e+05
eff_area[564]=8.0007981e+05&eff_area[565]=7.7814988e+05&eff_area[566]=7.5621994e+05
eff_area[567]=7.3429000e+05&eff_area[568]=7.1236006e+05&eff_area[569]=6.9043013e+05
eff_area[570]=6.6850019e+05&eff_area[571]=6.4686013e+05&eff_area[572]=6.2522007e+05
eff_area[573]=6.0358002e+05&eff_area[574]=5.8193996e+05&eff_area[575]=5.6029991e+05
eff_area[576]=5.3865985e+05&eff_area[577]=5.1701979e+05&eff_area[578]=4.9537974e+05
eff_area[579]=4.7373968e+05&eff_area[580]=4.5209962e+05&eff_area[581]=4.3263968e+05
eff_area[582]=4.1317973e+05&eff_area[583]=3.9371978e+05&eff_area[584]=3.7425984e+05
eff_area[585]=3.5479989e+05&eff_area[586]=3.3533994e+05&eff_area[587]=3.1588000e+05
eff_area[588]=2.9642005e+05&eff_area[589]=2.7696010e+05&eff_area[590]=2.5750016e+05
eff_area[591]=2.4143012e+05&eff_area[592]=2.2536008e+05&eff_area[593]=2.0929004e+05
eff_area[594]=1.9322001e+05&eff_area[595]=1.7714997e+05&eff_area[596]=1.6107993e+05
eff_area[597]=1.4500989e+05&eff_area[598]=1.2893986e+05&eff_area[599]=1.1286982e+05
eff_area[600]=9.6799781e+04&eff_area[601]=8.7173800e+04&eff_area[602]=7.7547819e+04
eff_area[603]=6.7921837e+04&eff_area[604]=5.8295856e+04&eff_area[605]=4.8669875e+04
eff_area[606]=3.9043894e+04&eff_area[607]=2.9417912e+04&eff_area[608]=1.9791931e+04
eff_area[609]=1.0165950e+04&eff_area[610]=5.3996863e+02&eff_area[611]=4.8597177e+02
eff_area[612]=4.3197490e+02&eff_area[613]=3.7797804e+02&eff_area[614]=3.2398118e+02
eff_area[615]=2.6998431e+02&eff_area[616]=2.1598745e+02&eff_area[617]=1.6199059e+02
eff_area[618]=1.0799373e+02&eff_area[619]=5.3996863e+01&eff_area[620]=0.0000000e+00
eff_area[621:899] = 0.0

end

; ---------------------------------------------------------------------
;
; End effect area subroutine 
;
; ---------------------------------------------------------------------



; ---------------------------------------------------------------------
;
; MAIN PROCEDURE
;
; ---------------------------------------------------------------------


pro dsc_advanced_kp_fit, year, doy, version, show=show, hold=hold, $
                         save=save, rezero=rezero, $
                         averaging_length = averaging_length, $
                         flybacks=flybacks, $
                         neg_offset=neg_offset, verbose=verbose, $
                         moment_threshold = moment_threshold,    $
                         clobber = clobber

mp = 1.67262178d-27             ; SI units
kb = 1.3806488d-23              ;  SI units
                       
; (1) INITIALIZATIONS AND LOADS
;compile necessary functions (Prchlik J. 2017/02/27)
;fist compile the startup program
;pref_set,'IDL_STARTUP','/crater/utilities/idl/mike/idlstartup',/commit
pref_set,'IDL_PATH','+/crater/utilities/idl/mike/:+/usr/local/itt/user_contrib/astron:+/usr/local/itt/user_contrib/coyote:+/usr/local/itt/user_contrib/catalyst:<IDL_DEFAULT>',/commit
load_plasma

;Grab needed cdf routines
setenv,'LD_LIBRARY_PATH=/usr/local/rsi/idl_6.3/bin/bin.linux.x86:/home/cdaweb/lib'
resolve_routine,'idlmakecdf',/COMPILE_FULL

;check if clobber is set
if keyword_set(clobber) then clober = 1 else clobber = 0

;Check to see if any file are created in current version (sav file created in idl_dscovr_to_cdf)
savefmt = '("dscovr_file_status_v",I02,".idl")'
savefil = string([version],format=savefmt)
version_save = file_test(savefil)

; (1.1) check if file is already process if version file exists when trying to save
if ((clobber ne 1) and keyword_set(save)) then begin



    ;Check on currect status of cdf file in version if version file exists
    if version_save eq 1 then begin
        ;restore the version savefile
        restore,savefil,/RELAXED_STRUCTURE_ASSIGNMENT
        ;check if file is fully processed (i.e. idl and cdf file already exists)
        cdf_find = size(where(((v_struct.obs_year eq year) and (v_struct.obs_doy eq doy))))
    ;if current version save file does not exits then will need to create a new save file so
    ;set cdf_find = 0
    endif else cdf_find=0

    ;Use idl file format to check if savefile exists 
    path =  '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute_corrected/'
    prefix =  '("dsc_fc_advkp_1minute_corrected_",I4,"_",I03,".idl")'
    idl_file = string([year,doy],format=prefix)
    idl_find = file_test(path+idl_file)

   
    

;Case to run only components if partially or fully processed
    case 1 of 
        idl_find eq 0: print,'Day Unprocessed Proceeding...'
        idl_find eq 1 and cdf_find[0] eq 0: begin
            print,'IDL file created but no cdf. Creating cdf'
            idl_dscovr_to_cdf,fix(year),fix(doy),version ;create 1minute cdf files based on doy and year
            compare_wind_dscovr,fix(year),fix(doy),version ;create plot comparing 1minute dscovr data to wind observations
            return
        end
        idl_find eq 1 and cdf_find[0] eq 1: begin
            print,'All files already created for given day (set clobber to overwrite).'
            return
        end
    endcase


endif



restore, '/crater/observatories/dscovr/code/modules/dsc_Vwindows.idl'
loadct, 39

; (1.2) Allows for xwindow plotting in batch mode
device, decomposed = 0
device, RETAIN=2
DLM_LOAD,'PNG' ;Added so PNG file creation cannot error

get_spectra, year, doy, ia, ib, ic, spec_jd, verbose=verbose ; get the current spectrum arrays
; these arrays ia, ib, and ic, are 64 x n where n is the number of
; time points in the day

jd = doy + julday(1, 1, year, 0, 0) - 1
ddd = string(doy, format = '(I03)')
yyyy = string(year, format = '(I4)')

if n_elements(moment_threshold) eq 0 then moment_threshold = 0.1

; (2) IDENTIFICATION OF CURRENT PEAKS
nspec = n_elements(spec_jd)
itot = ia+ib+ic
maxI = max(itot, dim=1, peaklocs)
hold_itot = itot
hold_ia = ia
hold_ib = ib
hold_ic = ic

; mark the step numbers and get speed guesses
nz = (hold_itot ne 0)
zeros = (hold_itot eq 0)
index = nz*((1.+fltarr(nspec)) ## indgen(64) + 1.)
fill = where(index eq 0)
index[fill] = 999.
firststep = min(index, dim = 1); this is the first step of the peak track cycle
trackstep = firststep+6 ; this is the MV_LO that the peak tracker keyed on
off_range = where(trackstep gt 63)
trackstep[off_range] = median(trackstep) ; fix any off range vals
trackstep = median(trackstep, 3) ; spike filter
v_guess = v_kms[trackstep]

index_array  =  (1.+fltarr(nspec)) ## findgen(64)
above_peak = index_array gt ((peaklocs mod 64) ## (1.+fltarr(64)))

; (3) REZEROING PROCEDURE
; re-zeroing procedure will use the minimum measured current from
; the high-energy side of the peak (i.e. avoid flyback)
;
; In the final version of the advkp, this may be made a default
; Make rezero the default state (Prchlik J. 2017/03/03)
if keyword_set(rezero) then rezero = rezero else rezero = 1

if n_elements(rezero) gt 0 then begin
  if rezero eq 1 then begin 
    if not keyword_set(neg_offset) then neg_offset = 0
    minIa = min(ia*above_peak + 999.*(ia*above_peak eq 0), dim = 1, minlocs_a)
    minIb = min(ib*above_peak + 999.*(ia*above_peak eq 0), dim = 1, minlocs_b)
    minIc = min(ic*above_peak + 999.*(ia*above_peak eq 0), dim = 1, minlocs_c)
    setpoints_a = median(ia[minlocs_a] < neg_offset, 3) ; spike filter and reject positive offsets
    setpoints_b = median(ib[minlocs_b] < neg_offset, 3) ; spike filter and reject positive offsets
    setpoints_c = median(ic[minlocs_c] < neg_offset, 3) ; spike filter and reject positive offsets
    offset_array_a = setpoints_a  ##  (1.+fltarr(64))
    offset_array_b = setpoints_b  ##  (1.+fltarr(64))
    offset_array_c = setpoints_c  ##  (1.+fltarr(64))
    ia = hold_ia - (hold_ia ne 0)*offset_array_a
    ib = hold_ib - (hold_ib ne 0)*offset_array_b
    ic = hold_ic - (hold_ic ne 0)*offset_array_c
    itot = ia+ib+ic  
 endif
endif

; This is our first derived set of MEASUREMENT uncertainties:
; -----------------------------------------------------
;|   offset_array_a, offset_array_b, offset_array_c    |
; -----------------------------------------------------


; if this is old data with flyback noise, zero out the flyback
; channels and proceed
if keyword_set(flybacks) then begin
  okstep = index_array  gt  ((firststep+2) ## (1.+fltarr(64)))
  fbstep = index_array  le  ((firststep+2) ## (1.+fltarr(64)))
  itot = itot*okstep
  ia=ia*okstep
  ib=ib*okstep
  ic=ic*okstep
  ia_FBERROR = ia*fbstep
  ib_FBERROR = ib*fbstep
  ic_FBERROR = ic*fbstep
endif else begin; if not, this uncertainty factor is just the minimum/most negative current
  ia_FBERROR = abs(min(ia, dim = 1)) ## (1.+fltarr(64)); a leading factor of 2 is a hack that
  ib_FBERROR = abs(min(ib, dim = 1)) ## (1.+fltarr(64)); produces uncertainties that are
  ic_FBERROR = abs(min(ic, dim = 1)) ## (1.+fltarr(64)); empirically consistent
endelse


; This is a second derived set of MEASUREMENT uncertainties
; but it usually won't be used.
; -----------------------------------------------------
;|   offset_array_a, offset_array_b, offset_array_c    |
;|   ia_FBERROR, ib_FBERROR, ic_FBERROR                |
; -----------------------------------------------------


; (3.5) calculate 1min averages. a 3pt median spike filter 
;       might also be a good idea somewhere in here, 
;       except we don't really see a lot of actual current spikes
if not keyword_set(averaging_length) then averaging_length = 15
ia1 = smooth(ia, [0, averaging_length])
ib1 = smooth(ib, [0, averaging_length])
ic1 = smooth(ic, [0, averaging_length])
itot1 = ia1+ib1+ic1

; (4) GET EFFECTIVE AREAS AND FLOW ANGLES
; (4.1) full resolution
maxI = max(itot, dim=1, peaklocs)
ABC2pt, ia[peaklocs], ib[peaklocs], ic[peaklocs], phi, theta, /corrected
alpha = sqrt(phi^2 + theta^2)     ; absolute angle with cup normal
aeff = dfc_geteffarea(phi, theta) ; effective area

; (4.2) 1min averages
maxI1 = max(itot1, dim=1, peaklocs1)
ABC2pt, ia1[peaklocs], ib1[peaklocs], ic1[peaklocs], phi1, theta1, /corrected
alpha1 = sqrt(phi1^2 + theta1^2)
aeff1 = dfc_geteffarea(phi1, theta1)

; (4.3) full res angular uncertainties
;       we're going to hack this from the offset correction
;       directly instead of calculating the curvature of the 
;       lookup table and propagating and all that. This should
;       be an upper bound.
ABC2pt, hold_ia[peaklocs], hold_ib[peaklocs], hold_ic[peaklocs], phis, thetas, /corrected
alphas = sqrt(phis^2 + thetas^2)

; uncertainty on phi and theta are based on comparison with the 
; pre-rezeroed calc
dphi = sqrt((phis - phi)^2)
dtheta = sqrt((thetas - theta)^2)
; uncertainty on the 1min averages include that and <1 min variance
; the latter is an uncertainty on the mean, so divide by sqrt(n)
phi1_sq = smooth(phi^2, averaging_length)
sigphi1 = sqrt(phi1_sq - smooth(phi, averaging_length)^2)/sqrt(averaging_length)
theta1_sq = smooth(theta^2, averaging_length)
sigtheta1 = sqrt(theta1_sq - smooth(theta, averaging_length)^2)/sqrt(averaging_length)
dphi1 = sqrt(dphi^2 + sigphi1^2)
dtheta1 = sqrt(dtheta^2 + sigtheta1^2)


; (5) CALCULATE THE FULL RES Reduced Distribution Functions
q0 = 1.6021892d-7 ; picocoulombs
eA = q0*aeff ## (1.+fltarr(64)); coulomb*cm^3/km
v = (1.+fltarr(nspec))##v_kms
dv = (1.+fltarr(nspec))##dv_kms
fa = ia/(eA*v*dv)
fb = ib/(eA*v*dv)
fc = ic/(eA*v*dv)
f = itot/(eA*v*dv) ; (picocoulomb * s-1)*km-2*s2 * km * cm-3 / picocoulomb
;                              = cm-3 * s * km-1 = 
;                              particles per cm-3 per km/s
; correct for NANs and infinities from bad i/o
tknan = where(finite(f, /nan) or finite(f, /inf), nnan)
f[tknan] = 0.
itot[tknan] = 0.

; (6) CALCULATE THE 1-minute average RDFs (is this a good approach, or
; should I average the currents?)
eA1 = q0*aeff1 ## (1.+fltarr(64))
f1 = itot1/(eA1*v*dv)

; (7) CALCULATE THE CHANNEL VARIANCES
m_fsq =  smooth(f^2, [0,averaging_length])
mf_sq = smooth(f, [0,averaging_length])^2
sigf1 = sqrt(m_fsq - mf_sq)


; Table of uncertainties tracked to this point
; -----------------------------------------------------------------------------
;|   ALL MEASUREMENTS:       offset_array_a, offset_array_b, offset_array_c   
;|   ALL MEASUREMENTS:       ia_FBERROR, ib_FBERROR, ic_FBERROR       
;|   1min RDF MEASUREMENTS:  sigf1
;|   flow angles:            dphi, dtheta
;|   1min flow angles:       dphi1, dtheta1
;|
;| other basic uncertainties
;| I had been using a base 5 pA uncertainty on all currents
;|    plus a f*0.02 ADC error
;|    plus a f*dv/v binning error
;|
; -----------------------------------------------------------------------------


; (8) Create variables for tracked uncertainties in the form df/f
df_base = 5./(eA*v*dv)
df_adc = f*0.02
df_binning = f*dv/v
df_offsets = sqrt(offset_array_a^2 + offset_array_b^2 + $
                        offset_array_c^2)/(eA*v*dv)
df_flyback = sqrt(ia_FBERROR^2 + ib_FBERROR^2 + $
                        ic_FBERROR^2)/(eA*v*dv)
df = sqrt(df_base^2 + df_adc^2 + df_binning^2 + df_flyback^2 + df_offsets^2)

df1_variance = sigf1
; treat offset, adc, and base uncertainties as un-correlated (may be a
; bad assumption for offset)
; The flybacks are generally quite correlated.
df1_adc = f1*0.02 /sqrt(averaging_length)
df1_offsets = smooth( sqrt(offset_array_a^2 + offset_array_b^2 + $
                        offset_array_c^2), [0, averaging_length])/ $
              (eA1*v*dv*sqrt(averaging_length))
df1_base = 5./(eA1*v*dv*sqrt(averaging_length))
df1_flyback = smooth(sqrt(ia_FBERROR^2 + ib_FBERROR^2 + $
                          ic_FBERROR^2), [0, averaging_length])/$
              (eA1*v*dv)
df1_binning = f1*(dv/v)/sqrt(averaging_length)
df1 = sqrt(df1_variance^2 + df1_adc^2 + df1_offsets^2 + $
           df1_base^2 + df1_flyback^2+df1_binning^2)

; (9) Create arrays for derived data
fitpars = fltarr(3, nspec)
fitsig = fltarr(3, nspec)
fitchisq = fltarr(nspec)
fit_nchan = intarr(nspec)
adv_fitpars = fltarr(6, nspec)
adv_chisq = fltarr(nspec)
adv_fitsig = fltarr(6, nspec)
adv_nchan = intarr(nspec)
moments = fltarr(4, nspec)

fitpars_1min = fltarr(3, nspec)
fitsig_1min = fltarr(3, nspec)
fitchisq_1min = fltarr(nspec)
fit_nchan_1min = intarr(nspec)
adv_fitpars_1min = fltarr(6, nspec)
adv_chisq_1min = fltarr(nspec)
adv_fitsig_1min = fltarr(6, nspec)
adv_nchan_1min = intarr(nspec)
moments_1min = fltarr(4, nspec)

; (10) Set up plot ranging
ok = where(f gt 0 and f lt 10, nok)
temp = f[ok]
stemp = sort(temp)
yrange = [-0.02, temp[stemp[0.99*nok]]]
xrange = [200, 1200]
xtickv = [200, 250, 300, 350, 400, 450, 500, 550, $
          600, 650, 700, 750, 800, 850, 900, 950, $
          1000, 1050, 1100, 1150, 1200]
xtickname =  ['200', ' ', '300', ' ','400',' ', '500', ' ','600', ' ','700', $
              ' ','800',' ', '900', ' ','1000', ' ','1100',' ', '1200']
xticks = n_elements(xtickv)-1

; Get spectrum timing in human units
caldat, spec_jd, thismonth, thisday, thisyear, hour, minute, second


; Initialize first loop and ascissa variables
i= 0L
hold_x0 = 30.
x = v_kms
dx = dv_kms

; (11) fitting loop --- run the fits on the 1min averages
;                       and then use the results for guesses
;                       on the full resolution


; This error handler skips the problematic spectrum and advances the
; loop index
CATCH, Error_status
                                ;This statement begins the error handler:
IF Error_status NE 0 THEN BEGIN
   PRINT, 'Error index: ', Error_status
   PRINT, 'Error message: ', !ERROR_STATE.MSG
   i=i+1
   PRINT, 'Skipping this spectrum'
   CATCH, /CANCEL
ENDIF

while i lt nspec do begin  
   
   ; (11.1) set up spectrum for fitting
   ; pull out this round of dep variables to fit
    y = reform(f1[*, i])  
    dy = reform(df1[*, i])
    ymax = max(y, x0); identify current peak channel
    w = 1./(dy^2) ; weight data points by 1/uncertainty^2 (std gaussian weight)

    ; (11.2) calculate moments
    momregions = label_region(y gt moment_threshold*ymax); select valid data region around peak
    momsel = where(y gt moment_threshold*ymax and momregions eq momregions[x0])
    if n_elements(momsel) le 1 then momsel = indgen(n_elements(x))
    mom0 = total((y*dx)[momsel])  
    mom1 = total((y*x*dx)[momsel])/mom0  
    mom2 = total((y*(x-mom1)*(x-mom1)*dx)[momsel])/mom0  
    mom3 = total((y*(x-mom1)*(x-mom1)*(x-mom1)*dx/mom0)[momsel])/mom0  
;    moments_1min[*, i] = [mom0, mom1, sqrt(abs(mom2)), (mom3/abs(mom3))*(abs(mom3))^(1./3.)]  
    ; emulate the NOAA algorithm moment correction, which 
    ; uses a maxwellian assumption to account for the truncated
    ; part of the distribution
    F_ratio = min(y[momsel])/max(y[momsel])
    qfac = sqrt(-alog(F_ratio))
    D0P = 1./erf(qfac)
    D1P = 1.
    D2P = 1./(erf(qfac) - 2*qfac*exp(-qfac^2)/sqrt(!dpi))
    nmom = d0p*mom0
    vmom = d1p*mom1
    tmom = d2p*1d6*(mp/kb)*mom2
    moments_1min[*,i] = [nmom, vmom, tmom, mom3]

    ; (11.3) calculate SIMPLE one gaussian fit
    if (x0 le 2) or (x0 ge 60) then x0 = hold_x0   
                                ; if this spec is messed up, try the
                                ; peak from last spec
    regions = label_region(y gt 0.025); select valid data region around peak
    sel = where(y gt 0.025 and regions eq regions[x0])
    if n_elements(sel) lt 5 then sel = [x0-2, x0-1, x0, x0+1, x0+2] 
    if (y[min(sel)-1] lt y[min(sel)]) and (y[min(sel)-1] gt 0) $
      then sel = [min(sel) - 1, sel]; left pad selection if y decreasing
    if (y[min(sel)-1] lt y[min(sel)]) and (y[min(sel)-1] gt 0) $
      then sel = [min(sel) - 1, sel]; left pad selection if y decreasing
    if max(sel) ge 63 then sel = sel-1   ; if running off range high translate selection
    ; use basic gaussfit for guess since this is usually robust
    fit = gaussfit(x[sel], y[sel], A, nterms = 3) 
    ; use weighted curve fit for improved 1-gaussian
    n_basic = n_elements(sel)
    fit = CALL_FUNCTION('curvefit', x[sel], y[sel], w[sel], A, sigma_basic, /DOUBLE, $ 
                        FUNCTION_NAME='one_gaussian', $
                       CHISQ=chisq_basic)     
    fitpars_1min[*, i] = a
    fitsig_1min[*, i] = sigma_basic
    fitchisq_1min[i] = chisq_basic
    fit_nchan_1min[i] = n_elements(sel)
 
    ; (11.4) calculate two-gaussian fit
    ; Construct a guess for the two-gaussian
    vv =  a[1]<v_guess[i] ; one of the two gaussians is always initialized
                          ; at or near the v_guess (the peak tracker value)
    ww = ((0.8*a[2]) > 5.) < 100.
    AA = [a[0], vv, ww, 0.4*a[0], vv+80., ww]  
    use = 1 - ((finite(y, /nan) or finite(y, /inf) or finite(dy, /nan) $
             or finite(y, /inf)))
    sel2 = where(use eq 1)

    test = CALL_FUNCTION('curvefit', x[sel2], y[sel2], w[sel2], aa, sigma, /DOUBLE, $ 
                        FUNCTION_NAME='two_gaussians', $
                       /NODERIVATIVE, ITMAX=itmax, ITER=it, tol=tol_init, CHISQ=chisq )  
    if abs(aa[5]) gt 250. then aa = [a[0], a[1], a[2], 0., 0., -9999.]  
    if abs(aa[2]) gt 250. then aa = [a[0], a[1], a[2], 0., 0., -9999.]  
    if aa[1] gt aa[4] and aa[4] ne 0. then begin  
       aaa = aa[0:2]  
       aa = [aa[3:5], aaa]  
    endif  

    adv_fitpars_1min[*, i] = aa
    adv_fitsig_1min[*, i] = sigma
    adv_chisq_1min[i] = chisq
    adv_nchan_1min[i] = n_elements(sel2)
   
    ; (11.5) optional real-time plotting of 1-minute fits
    if n_elements(show) gt 0 then begin   
      if show eq 1 then begin  
        plot, x, y, thick = 2, xrange = xrange, xstyle = 1, $
          psym = 10 ,  $
          ytitle ='AVERAGED RDF' , xticks = xticks, $
          xtitle = 'v [km/s]', ystyle = 1, /xl, xtickname=xtickname, $
          yrange = [yrange[0], yrange[1]>max(y)], xtickv=xtickv
        oploterr, x, y, dy, 3  
        oplot, x[sel], fit, thick = 5
        xx = 200 + 2.*findgen(500)  
        yy1 =  AA[0]*exp(-0.5* ((xx-AA[1])/AA[2])^2)  
        yy2 =  AA[3]*exp(-0.5* ((xx-AA[4])/AA[5])^2)  
        oplot, xx, yy1-0.0001, color= 254./4., thick = 2  
        oplot, xx, yy2-0.0001, color= 254./2., thick = 2  
        oplot, xx, yy1+yy2-0.0001, color = 254, thick = 3   
        oplot, x, -0.0001+0.04*2*AA[0]*exp(-0.5* ((x-sqrt(2.)*AA[1])/(0.5*AA[2]))^2) + $
           0.04*2*AA[3]*exp(-0.5* ((x-sqrt(2.)*AA[4])/(0.5*AA[5]))^2) , color = 3.*254./4., thick = 2  
        oplot, x[sel], y[sel], psym = 6, thick = 4, symsize = 0.5   
        xyouts, 10^!x.crange[1], !y.crange[1], align = 1, '!c!c' + string(hour[i], format = '(I02)') + $
                                                              ':' + string(minute[i], format = '(I02)') + $
                                                              ':' + string(second[i], format = '(I02)') + $
                                                              ' UT ', charsize = 1.   
        wait, 0.01  
     endif  
   endif; end of plotting conditional

    ; (11.6) now perform the same process on the high-resolution
   ; pull out this round of dep variables to fit
    y = reform(f[*, i])  
    dy = reform(df[*, i])
    ymax = max(y, x0); identify current peak channel
    w = 1./(dy^2) ; weight data points by 1/uncertainty^2 (std gaussian weight)

    ; (11.7) calculate moments
    momregions = label_region(y gt moment_threshold*ymax); select valid data region around peak
    momsel = where(y gt moment_threshold*ymax and momregions eq momregions[x0])
    ;if n_elements(momsel) le 1 then momsel = indgen(n_elements(x))
    if n_elements(momsel) lt 5 then momsel = [x0-2, x0-1, x0, x0+1, x0+2] ; mls 4/6 noaa pipeline selection (1)
    ok = where(momsel lt 63 and momsel gt 1) ; MLS 4/6 noaa pipeline selection (2)
    momsel = momsel[ok]                      ; MLS 4/6 noaa pipeline selection (3)

    yg0 = y>0 ; MLS 4/18 possible to have negative y on low energy side of the peak. Coordinated with Jeff
    mom0 = total((yg0*dx)[momsel])  
    mom1 = total((yg0*x*dx)[momsel])/mom0  
    mom2 = total((yg0*(x-mom1)*(x-mom1)*dx)[momsel])/mom0  
    mom3 = total((yg0*(x-mom1)*(x-mom1)*(x-mom1)*dx/mom0)[momsel])/mom0  
    F_ratio = min(yg0[momsel])/max(yg0[momsel])
    qfac = sqrt(-alog(F_ratio))
    D0P = 1./erf(qfac)
    D1P = 1.
    D2P = 1./(erf(qfac) - 2*qfac*exp(-qfac^2)/sqrt(!dpi))
    nmom = d0p*mom0
    vmom = d1p*mom1
    tmom = d2p*1d6*(mp/kb)*mom2
    moments[*,i] = [nmom, vmom, tmom, mom3]

    ; (11.8) calculate SIMPLE one gaussian fit
    if (x0 le 2) or (x0 ge 60) then x0 = hold_x0   
                                ; if this spec is messed up, try the
                                ; peak from last spec
    regions = label_region(y gt 0.025); select valid data region around peak
    sel = where(y gt 0.025 and regions eq regions[x0])
    if n_elements(sel) lt 5 then sel = [x0-2, x0-1, x0, x0+1, x0+2] 
    if (y[min(sel)-1] lt y[min(sel)]) and (y[min(sel)-1] gt 0) $
      then sel = [min(sel) - 1, sel]; left pad selection if y decreasing
    if (y[min(sel)-1] lt y[min(sel)]) and (y[min(sel)-1] gt 0) $
      then sel = [min(sel) - 1, sel]; left pad selection if y decreasing
    if max(sel) ge 63 then sel = sel-1   ; if running off range high translate selection
    ; use basic gaussfit for guess since this is usually robust
    fit = gaussfit(x[sel], y[sel], A, nterms = 3) 
    ; use weighted curve fit for improved 1-gaussian
    n_basic = n_elements(sel)
    fit = CALL_FUNCTION('curvefit', x[sel], y[sel], w[sel], A, sigma_basic, /DOUBLE, $ 
                        FUNCTION_NAME='one_gaussian', $
			CHISQ=chisq_basic)     
    fitpars[*, i] = a
    fitsig[*, i] = sigma_basic
    fitchisq[i] = chisq_basic
    fit_nchan[i] = n_elements(sel)
 
    ; (11.9) calculate two-gaussian fit
    ; Construct a guess for the two-gaussian
    vv =  a[1]<v_guess[i] ; one of the two gaussians is always initialized
                          ; at or near the v_guess (the peak tracker value)
    ww = ((0.8*a[2]) > 5.) < 100.
    AA = [a[0], vv, ww, 0.4*a[0], vv+80., ww]  
    use = 1 - ((finite(y, /nan) or finite(y, /inf) or finite(dy, /nan) $
             or finite(y, /inf)))
    sel2 = where(use eq 1)

    test = CALL_FUNCTION('curvefit', x[sel2], y[sel2], w[sel2], aa, sigma, /DOUBLE, $ 
                        FUNCTION_NAME='two_gaussians', $
			/NODERIVATIVE, ITMAX=itmax, ITER=it, tol=tol_init, CHISQ=chisq )  
    if abs(aa[5]) gt 250. then aa = [a[0], a[1], a[2], 0., 0., -9999.]  
    if abs(aa[2]) gt 250. then aa = [a[0], a[1], a[2], 0., 0., -9999.]  
    if aa[1] gt aa[4] and aa[4] ne 0. then begin  
       aaa = aa[0:2]  
       aa = [aa[3:5], aaa]  
    endif  

    adv_fitpars[*, i] = aa
    adv_fitsig[*, i] = sigma
    adv_chisq[i] = chisq
    adv_nchan[i] = n_elements(sel2)
    ; (11.10) optional real-time plotting of high-res fits
    if n_elements(show) gt 0 then begin   
      if show ne 1 then begin  
        plot, x, y, thick = 2, xrange = xrange, xstyle = 1, $
          psym = 10 ,  $
          ytitle ='RDF' , xticks = xticks, $
          xtitle = 'v [km/s]', ystyle = 1, /xl, xtickname=xtickname, $
          yrange = [yrange[0], yrange[1]>max(y)], xtickv=xtickv
        oploterr, x, y, dy, 3  
        oplot, x[sel], fit, thick = 5
        xx = 200 + 2.*findgen(500)  
        yy1 =  AA[0]*exp(-0.5* ((xx-AA[1])/AA[2])^2)  
        yy2 =  AA[3]*exp(-0.5* ((xx-AA[4])/AA[5])^2)  
        oplot, xx, yy1-0.0001, color= 254./4., thick = 2  
        oplot, xx, yy2-0.0001, color= 254./2., thick = 2  
        oplot, xx, yy1+yy2-0.0001, color = 254, thick = 3   
        oplot, x, -0.0001+0.04*2*AA[0]*exp(-0.5* ((x-sqrt(2.)*AA[1])/(0.5*AA[2]))^2) + $
           0.04*2*AA[3]*exp(-0.5* ((x-sqrt(2.)*AA[4])/(0.5*AA[5]))^2) , color = 3.*254./4., thick = 2  
        oplot, x[sel], y[sel], psym = 6, thick = 4, symsize = 0.5   
        xyouts, 10^!x.crange[1], !y.crange[1], align = 1, '!c!c' + string(hour[i], format = '(I02)') + $
                                                              ':' + string(minute[i], format = '(I02)') + $
                                                              ':' + string(second[i], format = '(I02)') + $
                                                              ' UT ', charsize = 1.   
        wait, 0.01  
     endif  
   endif; end of plotting conditional
 
    ; prep for next iteration
    hold_x0 = x0  
    i=i+1
 endwhile ; END OF FITTING LOOP


; (12) store all of the fitting info into a structure for saving
adv_kp =  {fitpars:fitpars, $
           fitchisq:fitchisq, $
           fitsig:fitsig, $
           fit_nchan:fit_nchan, $
           adv_fitpars:adv_fitpars, $
           adv_chisq:adv_chisq, $
           adv_fitsig:adv_fitsig, $
           adv_nchan:adv_nchan, $
           fitpars_1min:fitpars_1min, $
           fitchisq_1min:fitchisq_1min, $
           fitsig_1min:fitsig_1min, $
           fit_nchan_1min:fit_nchan_1min, $
           adv_fitpars_1min:adv_fitpars_1min, $
           adv_chisq_1min:adv_chisq_1min, $
           adv_fitsig_1min:adv_fitsig_1min, $
           adv_nchan_1min:adv_nchan_1min, $
           moments:moments, $
           moments_1min:moments_1min, $
           phi:phi, dphi:dphi, $
           theta:theta, dtheta:dtheta, $
           phi_1min:phi1, dphi_1min:dphi1, $
           theta_1min:theta1, dtheta_1min:dtheta1, $
           spec_jd:spec_jd}

; (13) translate highres fitting results into human-readable key parameters

; BASIC FITPAR EXTRACTION AND PREP: Full resolution
notes =  ['Prepared with dsc_advanced_kp_fit_v2', $
                   'FULL RESOLUTION nonlinear fits', $
                   'Single maxwellian model', $
                   'active background subtraction', $
                   'Post-optimization data. IT=21', $
                   'science quality.', $
                   'mstevens@cfa.harvard.edu '+systime()] 
fit_output_structure, fitpars, fitsig, phi, theta, dphi, dtheta, year, spec_jd, dfc_kp, notes = notes
rangechecks, dfc_kp

; 1 minute average fitpar extraction and prep
notes =  ['Prepared with dsc_advanced_kp_fit_v2', $
                   'nonlinear fits to 1-minute averaged spectra', $
                   'Single maxwellian model', $
                   'active background subtraction', $
                   'Post-optimization data. IT=21', $
                   'science quality.', $
                   'mstevens@cfa.harvard.edu '+systime()]

; use regridding routine to downsample the 1 minute data
invars = [[transpose(fitpars_1min)], $
          [transpose(fitsig_1min)], $
          [phi1], $
          [theta1], $
          [dphi1], $
          [dtheta1]]
thisdoy = spec_jd - julday(1, 1, year, 0, 0, 0) + 1.
regrid_1minute, invars, outvars, doy, thisdoy, outdoy=outdoy
gridjd = outdoy + julday(1, 1, year, 0, 0, 0) -1
fit_output_structure, transpose(outvars[*, 0:2]), transpose(outvars[*, 3:5]), outvars[*, 6], $
                      outvars[*, 7], outvars[*, 8], outvars[*, 9], year, gridjd, dfc_kp_1min, notes = notes
rangechecks, dfc_kp_1min
                            
; (14) Save if requested, hault execution if requested
if keyword_set(save) then begin ; (Prchlik. J 2017/02/27 added cdf output and compare plot under save)
    advkp_save, adv_kp, dfc_kp, dfc_kp_1min, yyyy, ddd
    ;load newly created save file and apply correction
    apply_empirical_corrections_1min, fix(yyyy),fix(ddd)-45,fix(ddd)+45
    idl_dscovr_to_cdf,fix(yyyy),fix(ddd),version ;create 1minute cdf files based on doy and year
    compare_wind_dscovr,fix(yyyy),fix(ddd),version ;create plot comparing 1minute dscovr data to wind observations
;setup idl structures




endif

if keyword_set(hold) then stop


end


; ----------------------- WRAPPERS --------------------------------------




; ---------------------------------------------------------------------
;
; Procedure to process a range of days
;
; ---------------------------------------------------------------------

pro advkp_doyrange, year, d1, d2, version, save=save, rezero=rezero, show=show, clobber=clobber

for thisday = d1, d2 do dsc_advanced_kp_fit, year, thisday, version, save=save, rezero=rezero, show=show, clobber=clobber

end


; ---------------------------------------------------------------------
;
; Procedure to load 1-minute data for a range of days
; by keyword
;
; ---------------------------------------------------------------------

pro load_advkp_1min, year, d1, d2, doy=doy, n=n, vx=vx, vy=vy, vz=vz, t=t, w=w, $
                     dn=dn, dvx=dvx, dvy=dvy, dvz=dvz, dt=dt, dw=dw, corrected=corrected

yyyy = string(year, format = '(I4)')
prefix =  'dsc_fc_advkp_1minute'
if not keyword_set(corrected) then begin
   path = '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute/'
   prefix = 'dsc_fc_advkp_1minute'
endif else begin
   path =  '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute_corrected/'
   prefix =  'dsc_fc_advkp_1minute_corrected'
endelse

datfiles = file_search(path+ prefix+'_' + yyyy + '*.idl')
      
; exit if there are no data files found
if n_elements(datfiles) eq 1 and datfiles[0] eq '' then begin 
   path = dialog_pickfile(title = 'Please specify local 1 minute advkp directory', /dir)
   datfiles = file_search(new_dir+prefix+'_'+yyyy+'*.idl')
   if n_elements(datfiles) eq 1 and attfiles[0] eq '' then begin
      print, 'No data files found. load_advkp_1min failed.'
      return
   endif
endif

fileddd = strmid(datfiles, 6, 3, /reverse)
filedoy = fix(fileddd)
tk = where(filedoy ge (d1 - 1) and filedoy le (d2 + 1), ntk)
; exit if there are no attitude files in this range
if ntk eq 0 then begin 
   print, 'load_advkp_1min: no data files found in requested range.'
   return
endif

tkfiles = datfiles[tk]
restore, tkfiles[0]

vx = DFC_KP_1MIN.vx.data
vy = DFC_KP_1MIN.vy.data
vz = DFC_KP_1MIN.vz.data
dvx = DFC_KP_1MIN.vx.uncertainty
dvy = DFC_KP_1MIN.vy.uncertainty
dvz = DFC_KP_1MIN.vz.uncertainty
n = DFC_KP_1MIN.n.data
w = DFC_KP_1MIN.w.data
T = DFC_KP_1MIN.T.data
dn = DFC_KP_1MIN.n.uncertainty
dw = DFC_KP_1MIN.w.uncertainty
dT = DFC_KP_1MIN.T.uncertainty
doy = DFC_KP_1MIN.time.data

for i = 1, ntk-1 do begin &$
   restore, tkfiles[i] &$
   vx = [ vx, DFC_KP_1MIN.vx.data]
   vy = [ vy, DFC_KP_1MIN.vy.data]
   vz = [ vz, DFC_KP_1MIN.vz.data]
   dvx = [ dvx, DFC_KP_1MIN.vx.uncertainty]
   dvy = [ dvy, DFC_KP_1MIN.vy.uncertainty]
   dvz = [ dvz, DFC_KP_1MIN.vz.uncertainty]
   n = [ n, DFC_KP_1MIN.n.data]
   w = [ w, DFC_KP_1MIN.w.data]
   T = [ T, DFC_KP_1MIN.T.data]
   dn = [ dn, DFC_KP_1MIN.n.uncertainty]
   dw = [ dw, DFC_KP_1MIN.w.uncertainty]
   dT = [ dT, DFC_KP_1MIN.T.uncertainty]
   doy = [ doy, DFC_KP_1MIN.time.data]
endfor

tk = where(doy ge d1 and doy le d2, ntk)
if ntk gt 0 then begin
   vx = vx[tk]
   vy = vy[tk]
   vz = vz[tk]
   dvx = dvx[tk]
   dvy = dvy[tk]
   dvz = dvz[tk]
   n = n[tk]
   dn = dn[tk]
   w = w[tk]
   dw = dw[tk]
   T = T[tk]
   dT = dT[tk]
   doy = doy[tk]
endif

end



; ---------------------------------------------------------------------
;
; Procedure to derive a density correction factor/function
; using Wind SWE
;
; USES: load_plasma.pro
; @'/crater/utilities/idl/mike/idlstartup'
; ---------------------------------------------------------------------

pro derive_empirical_corrections, year, d1, d2, density_polynomial, wfactor

load_advkp_1min, year, d1, d2, doy=doy, n=n, vx=vx, vy=vy, vz=vz, t=t, w=w, $
                     dn=dn, dvx=dvx, dvy=dvy, dvz=dvz, dt=dt, dw=dw

load_plasma, year, d1, d2, /wind, doy=doyw, vx = vxw, vy=vyw, vz=vzw, den = nw, t= ww


; interpolates
u = -vx
uw = -vxw
uw_i = interpol(median(uw, 3), doyw, doy)
nw_i = interpol(median(nw, 3), doyw, doy)
ww_i = interpol(median(ww, 3), doyw, doy)

match = where(abs(uw_i - u) lt 10 and nw_i gt 0 and n gt 0 $
              and nw_i gt 0 and ww_i gt 0 and w gt 0 and n lt 100 and $
              (ww_i/w) gt 0.5 and (ww_i/w) lt 2., nmatch)

 dev = stddev( nw_i[match]/n[match])
 vr = u[match]                                            
 ratio = nw_i[match]/n[match]      
 ord = sort(vr)                        
 vr = vr[ord]
 ratio = ratio[ord]
 ratio = smooth(ratio, 60)
 vr = [vr, 900, 1000, 1100, 1200]
 ratio = [ratio, 1, 1, 1, 1]
 coeffs = poly_fit(vr, ratio, 6, yfit = yfit)
 plot, vr, smooth(ratio, 60), ytitle = 'Density ratio (Wind/DSC)', xtitle = 'DSC speed'
 oplot, vr, yfit, color = 254., thick =4

 density_polynomial = coeffs

; rescale and examine the dscovr densities
 yfac = (u^0.) * coeffs[0]
 for i = 1, n_elements(coeffs)-1 do  yfac = yfac + coeffs[i]*u^i 
 yfac = yfac*(u le 750.)*(u ge 250.) + (u lt 250) + (u gt 750)
 n_rescaled = n*yfac

; the dscovr thermal speed may be consistently ~5% low
 match = where(abs(uw_i - u) lt 10 and nw_i gt 0 and n gt 0)
 wfactor = median( (ww_i/w)[match] )



end



; ---------------------------------------------------------------------
;
; Procedure to apply empirical corrections from Wind SWE
;
; one we may still need is to put absolute limits on temperature,
; since there are still some spurious features getting fit every once
; in a while
;
; ---------------------------------------------------------------------

pro apply_empirical_corrections_1min, year, d1, d2, path=path, verbose=verbose,derive=derive
  
  if keyword_set(verbose) then print, 'deriving empirical corrections'
  derive_empirical_corrections, year, d1, d2, coeffs, wfactor 
  
  yyyy = string(year, format = '(I4)')
  prefix =  'dsc_fc_advkp_1minute'
  if not keyword_set(path) then path = '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute/'
  datfiles = file_search(path+ prefix+'_' + yyyy + '*.idl')
  outpath =  '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute_corrected/'
  outprefix =  'dsc_fc_advkp_1minute_corrected'
  
; exit if there are no data files found
  if n_elements(datfiles) eq 1 and datfiles[0] eq '' then begin 
     path = dialog_pickfile(title = 'Please specify local 1 minute advkp directory', /dir)
     datfiles = file_search(new_dir+prefix+'_'+yyyy+'*.idl')
     if n_elements(datfiles) eq 1 and attfiles[0] eq '' then begin
        print, 'No data files found. load_advkp_1min failed.'
        return
     endif
  endif
  
  fileddd = strmid(datfiles, 6, 3, /reverse)
  filedoy = fix(fileddd)
  tk = where(filedoy ge (d1 - 1) and filedoy le (d2 + 1), ntk)
; exit if there are no attitude files in this range
  if ntk eq 0 then begin 
     print, 'load_advkp_1min: no data files found in requested range.'
     return
  endif
  
  tkfiles = datfiles[tk]
  tkfiles_ddd = fileddd[tk]
  
  for thisfile = 0, ntk-1 do begin
     restore, tkfiles[thisfile]
; rescale and examine the dscovr densities
     u = -DFC_KP_1MIN.vx.data
     yfac = (u^0.) * coeffs[0]
     for i = 1, n_elements(coeffs)-1 do  yfac = yfac + coeffs[i]*u^i 
                                ; this next step is just a mask to
                                ; make sure only in-range, non-fill
                                ; data values are rescaled
     valid = where((DFC_KP_1min.w.data ne DFC_KP_1min.w.fillval) and $
                   (DFC_KP_1min.T.data ne DFC_KP_1min.T.fillval) and $
                   (DFC_KP_1min.n.data ne DFC_KP_1min.n.fillval) and $
                   (u ge 250.) and (u le 750), nvalid)
     if nvalid gt 0 then begin
        DFC_KP_1min.n.data[valid] = DFC_KP_1min.n.data[valid]*yfac[valid]
        DFC_KP_1min.w.data[valid] = DFC_KP_1min.w.data[valid]*wfactor
        DFC_KP_1min.T.data[valid] = DFC_KP_1min.T.data[valid]*wfactor^2
     endif
     temp = DFC_KP_1min.notes[1]
     DFC_KP_1min.notes[1] = temp + ', empirically rescaled with Wind SWE'
     ; re-apply range checks
     rangechecks, dfc_kp_1min
     outfile =  outpath + outprefix + '_' + yyyy + '_' + tkfiles_ddd[thisfile] + '.idl'
     save, DFC_KP_1min, filename = outfile
     if keyword_set(verbose) then $
        print, 'corrected ' +  outprefix + '_' + yyyy + '_' + tkfiles_ddd[thisfile] + '.idl'
  endfor


end

pro apply_empirical_corrections_fullres, year, d1, d2, path=path, verbose=verbose
  
  if keyword_set(verbose) then print, 'deriving empirical corrections'
  derive_empirical_corrections, year, d1, d2, coeffs, wfactor
  
  yyyy = string(year, format = '(I4)')
  prefix =  'dsc_fc_advkp_fullres'
  if not keyword_set(path) then path = '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/fullres/'
  datfiles = file_search(path+ prefix+'_' + yyyy + '*.idl')
  outpath =  '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/fullres_corrected/'
  outprefix =  'dsc_fc_advkp_fullres_corrected'
  
; exit if there are no data files found
  if n_elements(datfiles) eq 1 and datfiles[0] eq '' then begin 
     path = dialog_pickfile(title = 'Please specify local fullres advkp directory', /dir)
     datfiles = file_search(new_dir+prefix+'_'+yyyy+'*.idl')
     if n_elements(datfiles) eq 1 and attfiles[0] eq '' then begin
        print, 'No data files found. load_advkp_fullres failed.'
        return
     endif
  endif
  
  fileddd = strmid(datfiles, 6, 3, /reverse)
  filedoy = fix(fileddd)
  tk = where(filedoy ge (d1 - 1) and filedoy le (d2 + 1), ntk)
; exit if there are no files in this range
  if ntk eq 0 then begin 
     print, 'load_advkp: no data files found in requested range.'
     return
  endif
  
  tkfiles = datfiles[tk]
  tkfiles_ddd = fileddd[tk]
  
  for thisfile = 0, ntk-1 do begin
     restore, tkfiles[thisfile]
; rescale and examine the dscovr densities
     u = -DFC_KP.vx.data
     yfac = (u^0.) * coeffs[0]
     for i = 1, n_elements(coeffs)-1 do  yfac = yfac + coeffs[i]*u^i 
                                ; this next step is just a mask to
                                ; make sure only in-range, non-fill
                                ; data values are rescaled
     valid = where((DFC_KP.w.data ne DFC_KP.w.fillval) and $
                   (DFC_KP.T.data ne DFC_KP.T.fillval) and $
                   (DFC_KP.n.data ne DFC_KP.n.fillval) and $
                   (u ge 250.) and (u le 750), nvalid)
     if nvalid gt 0 then begin
        DFC_KP.n.data[valid] = DFC_KP.n.data[valid]*yfac[valid]
        DFC_KP.w.data[valid] = DFC_KP.w.data[valid]*wfactor
        DFC_KP.T.data[valid] = DFC_KP.T.data[valid]*wfactor^2
     endif
     temp = DFC_KP.notes[1]
     DFC_KP.notes[1] = temp + ', empirically rescaled with Wind SWE'
     ; re-apply range checks
     rangechecks, dfc_kp
     outfile =  outpath + outprefix + '_' + yyyy + '_' + tkfiles_ddd[thisfile] + '.idl'
     save, DFC_KP, filename = outfile
     if keyword_set(verbose) then $
        print, 'corrected ' +  outprefix + '_' + yyyy + '_' + tkfiles_ddd[thisfile] + '.idl'
  endfor


end


; ---------------------------------------------------------------------
;
; Summary plot procedure
;
; ---------------------------------------------------------------------

pro advkp_1min_corrected_summary_plot, year, day
 
  yyyy = string(year, format = '(I4)')
  ddd = string(day, format = '(I03)')
  prefix =  'dsc_fc_advkp_1minute_corrected'
  if not keyword_set(path) then $
     path = '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute_corrected/'
  datfile = file_search(path+ prefix+'_' + yyyy + '_' + ddd + '.idl')
  if n_elements(datfile) eq 1 and datfile[0] eq '' then $
     datfile = dialog_pickfile(title = 'Please specify advkp data file')
  restore, datfile
  filename =  '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute_corrected_plots/'+ $
          prefix + '_plot_'+yyyy+'_'+ddd+'.eps'

  set_plot, 'ps'
  device, /color, /inches, /helvetica, /encapsulate, xsize = 8.5, ysize = 11, $
          filename = filename
  !p.font = 0
  dat= DFC_KP_1MIN
  valid = where(dat.n.data ne dat.n.fillval)
  yoff = max( [abs(dat.vz.data[valid])+dat.vz.uncertainty[valid], $
               abs(dat.vy.data[valid])+dat.vy.uncertainty[valid]]) < 200
  pwin = [0.1, 0, 0.95, 0.15]
  dy = [0, 0.01, 0, 0.01]
  py = [0, 0.15, 0, 0.15]

  plot, dat.time.data[valid], dat.n.data[valid], position = pwin+10*dy, $
        ytitle = dat.n.name + '!c['+dat.n.units+']', $
        xtitle = dat.time.units + ', '+yyyy, psym = 3
  oploterror, dat.time.data[valid],  dat.n.data[valid], $
              dat.n.uncertainty[valid], /nohat, errcol = 'grey', psym = 3
  oplot, dat.time.data[valid], dat.n.data[valid], psym = 6, symsize = 0.15
  
  plot, dat.time.data[valid], dat.w.data[valid], position = pwin+10*dy+py+dy, $
        ytitle = 'proton thermal speed' + '!c['+dat.w.units+']', $
        /noerase, xtickformat = '(A1)', psym = 3
  oploterror, dat.time.data[valid],  dat.w.data[valid], $
              dat.w.uncertainty[valid], /nohat, errcol = 'grey', psym = 3
  oplot, dat.time.data[valid], dat.w.data[valid], psym = 6, symsize = 0.15

  plot, dat.time.data[valid], dat.vz.data[valid], position = pwin+10*dy+2.*(py+dy), $
        ytitle = 'proton v!dz,GSE!n' + '!c['+dat.vz.units+']', $
        /noerase, xtickformat = '(A1)', yrange = [-yoff, yoff], psym = 3
  oplot, !x.crange, [0,0]
  oploterror, dat.time.data[valid],  dat.vz.data[valid], $
              dat.vz.uncertainty[valid], /nohat, errcol = 'grey', psym = 3
  oplot, dat.time.data[valid], dat.vz.data[valid], psym = 6, symsize = 0.15
 
  plot, dat.time.data[valid], dat.vy.data[valid], position = pwin+10*dy+3.*(py+dy), $
        ytitle = 'proton v!dy,GSE!n' + '!c['+dat.vz.units+']', $
        /noerase, xtickformat = '(A1)', yrange = [-yoff, yoff], psym = 3
  oplot, !x.crange, [0,0]
  oploterror, dat.time.data[valid],  dat.vy.data[valid], $
              dat.vy.uncertainty[valid], /nohat, errcol = 'grey', psym = 3
  oplot, dat.time.data[valid], dat.vy.data[valid], psym = 6, symsize = 0.15

  plot, dat.time.data[valid], -dat.vx.data[valid], position = pwin+10*dy+4.*(py+dy), $
        ytitle = 'proton -v!dx,GSE!n' + '!c['+dat.vx.units+']', $
        /noerase, xtickformat = '(A1)', yrange = minmax( -dat.vx.data[valid]), psym = 3
  oploterror, dat.time.data[valid],  -dat.vx.data[valid], $
              dat.vx.uncertainty[valid], /nohat, errcol = 'grey', psym = 3
  oplot, dat.time.data[valid], -dat.vx.data[valid], psym = 6, symsize = 0.15

device, /close
set_plot, 'x'

end

   
; --------------------------- BEGIN SANDBOX ------------------------------


pro advkp_1min_corrected_makeCDF, year, doy

; identify idl file for CDF
  yyyy = string(year, format = '(I4)')
  ddd = string(day, format = '(I03)')
  prefix =  'dsc_fc_advkp_1minute_corrected'
  if not keyword_set(path) then $
     path = '/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute_corrected/'
  datfile = file_search(path+ prefix+'_' + yyyy + '_' + ddd + '.idl')
  if n_elements(datfile) eq 1 and datfile[0] eq '' then $
     datfile = dialog_pickfile(title = 'Please specify advkp data file')
  restore, datfile
  
; identify ISTP CDF skeleton file
;;; ... left this off on 1/30 to address flow angle issues. Will pick
;;; back up

end

pro scratch
end


;Allows you to compile the full file
pro dsc_advanced_kp_fit_v2

end
