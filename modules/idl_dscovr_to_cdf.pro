;idl_dscovr_to_cdf
;
;--------------------------------------------------
;Uses
;    IF not ran from dsc_advanced_kp_fit_v2 then the follow contingencies are required
;    following 
;    @'/crater/utilities/idl/mike/idlstartup' 
;    @compile_cdaweb
;    @compile_IDLmakecdf

;--------------------------------------------------
;ROUTINE for coverting idlsave file to CDF file
;Coverts idl save files to CDF format
;
;USAGE
;idl_dscovr_to_cdf,year,doy,version,archive=archive,skeleton=skeleton,filefmt=filefmt, $
;                  outfmt=outfmt,orchive=orchive,outdom=outdom
;Added keyword to output cdf in dom format
;--------------------------------------------------

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
;ldoy = leap_check(year)
;
;COMMENT
;Checks for leap year
;--------------------------------------------------
function leap_check,year
;leap year check
four = year/4.
;last day of year
ldoy = 365
if four-fix(four) lt .20 then ldoy = 366
return,ldoy
end

;--------------------------------------------------
;USAGE
;dtime = chi_min_time(wdoy,wd,ddoy,dd,ddd,span=span,samp=samp)
;
;COMMENTS
;    Function which spans +/- 10 minutes to find Chi^2 min
;    wdoy = WIND day of year
;    wd   = WIND-density 
;    ddoy = DSCOVR day of year
;    dd   = DSCOVR-density
;    ddd  = DSCOVR uncertainty in density
;    span = time in minutes to minimize
;    samp = sampling fequence in minutes
;--------------------------------------------------

function chi_min_time,wdoy,wd,ddoy,dd,ddd,span=span,samp=samp
if keyword_set(span) then span = span else span = 15 ;span to loop +/- in minutes
if keyword_set(samp) then samp = samp else samp = 0.50 ; sampling fequency in minutes

span = span/60./24. ;turn into fraction of a day
samp = samp/60./24. ;turn into fraction of a day

;array of offset to loop over
nbins = round(2.*span/samp)+1
tvals = findgen(nbins)*samp-span
cvals = fltarr(nbins) ; array of zeros to store chi^2 vals

;loop over points for chi2 sampling
for i=0,nbins-1 do begin
    ;create interpolation from wind data assuming small sigma
    sdd = interpol(wd,wdoy,ddoy+tvals[i])-dd
    chi2 = total(sdd^2/ddd^2)
;store chi2 val
    cvals[i] = chi2
endfor

;fine sample chi2 distirbution to find minimum value
fsamp = samp/100.
;array of offset to loop over
fnbins = round(2.*span/fsamp)+1
ftvals = findgen(fnbins)*fsamp-span
fcvals = interpol(cvals,tvals,ftvals,/LSQUADRATIC) ; array of of chi^2 vals

;minimum chi^2 
chim = where(fcvals eq min(fcvals))
if n_elements(size(chim)) gt 3 then begin ;move the time as little as possible
    dtime = ftvals[chim]
    minmo = where(abs(dtime) eq min(abs(dtime))) ; move the minimum amount
    dtime = dtime[minmo]
endif else begin 
    dtime = 0. ; return 0 if no min??
endelse

print,'Chi^2 min value',strcompress(min(fcvals),/remove_all)
print,'Chi^2 time offset',strcompress(dtime*24.*60.),'min'


return,dtime
end

;--------------------------------------------------
;
;USAGE
;grp = grp_obs(y,npix)
;
;COMMENTS
;Returns an array of indices which demark the beginning and end of a group
;--------------------------------------------------
function grp_obs,y,npix

ind_grp = [0]

;looks at each interior element if the diff is 1 contiue else add y[i] and y[i+1] do array
for i=1,n_elements(y)-2 do if y[i+1]-y[i] ne 1 then ind_grp = [temporary(ind_grp), i, i+1]

;include the last element always
ind_grp = [temporary(ind_grp),n_elements(y)-1]

g1 = ind_grp[0:n_elements(ind_grp)-2:2]
g2 = ind_grp[1:n_elements(ind_grp)-1:2]

;get where the grouped array is greater than the npix values
bad = where(g2-g1 ge npix,bcnt)

;if regions where npix different from npix exist compute times
if bcnt gt 0 then begin
    ;starting and ending bad regions in ind_grp array
    s1 = g1[bad]
    e1 = g2[bad]
    
    ;starting and ending bad regions in y array
    s2 = y[s1]
    e2 = y[e1]
    
    
    ;covert split arrays into indices array
    ;also start the bad pixel npix before and end npix after
    ;for i=0,n_elements(s2)-1 do if i eq 0 then grp = findgen(fix(e2[i]-s2[i])+1)+s2[i] else grp = [temporary(grp),findgen(fix(e2[i]-s2[i])+1)+s2[i]]
    for i=0,n_elements(s2)-1 do begin
        s_cut = s2[i]-npix 
        e_cut = e2[i]+npix 
        ;prevent s_cut or e_cut from going out of bounds
        if s_cut lt 0 then s_cut=0
        if e_cut gt 1440 then e_cut=1440
        print,s_cut,e_cut
        if i eq 0 then grp = findgen(fix(e_cut-s_cut)+1)+s_cut else grp = [temporary(grp),findgen(fix(e_cut-s_cut)+1)+s_cut]
    endfor
endif else grp = -9999.0


return,grp
end




;--------------------------------------------------
;
;USAGE
;replace = den_flag(den_val,den_err,doy_val,sigcut=sigcut,npix=npix)
;
;COMMENTS
;Sends flag for densities which differ by sigcut*measures_sigma from WIND
;
;--------------------------------------------------
function den_flag,den_val,den_unc,doy_val,year,sigcut=sigcut,npix=npix
if keyword_set(npix) then npix = npix else npix = 25 ; 25 minute check
if keyword_set(sigcut) then sigcut=sigcut else sigcut = 5.00


;get doy 
doy = fix(doy_val[0])

;start doy
sdoy = doy - 1.
;end doy
edoy = doy + 2.

ldoy = leap_check(year)


case 1 of
    (sdoy lt 1): begin
;If first doy start search with last day of previous year
        syear = year-1
        sdoy = leap_check(year)
        load_plasma,syear,sdoy,sdoy+1,/wind,doy=wdoy1,den=wden1
        load_plasma,year,1,edoy,/wind,doy=wdoy2,den=wden2
        wdoy = [wdoy1-sdoy,wdoy2]
        wden = [wden1,wden2]
    end
;If last doy end search with first day of next year        
    (edoy gt ldoy+1): begin
        load_plasma,year,doy,sdoy,/wind,doy=wdoy1,den=wden1
        load_plasma,year+1,1,2,/wind,doy=wdoy2,den=wden2
        wdoy = [wdoy1,wdoy2+sdoy]
        wden = [wden1,wden2]
    end
;No year break just use load plasma
    else: load_plasma,year,sdoy,edoy,/wind,doy=wdoy,den=wden
endcase
     
;time offset between WIND and DSCOVR
dtime = chi_min_time(wdoy,wden,doy_val,den_val,den_unc,span=span,samp=samp)

;Difference between observed DSCOVR and interpolated WIND data
wind_den = interpol(wden,wdoy,doy_val+dtime[0])
chk_den = where(den_val gt -9990.)
del_den = (den_val[chk_den]-wind_den[chk_den])/wind_den[chk_den]
;sig_off = den_unc[chk_den]/den_val[chk_den]

;Calculate the moving average offset and measurement error
del_off = ts_smooth(del_den,npix)
sig_emp = ts_sigma(del_den,npix) ;emperical scatter
sig_off = ts_smooth(den_unc[chk_den]/den_val[chk_den],npix)

;combined measurement error and emperical uncertainty
sig_off = sqrt(sig_emp^2+temporary(sig_off)^2)/sqrt(2.)

;Do a 3 sigma rejection areas where average offset and uncertainty does not cross 0
bad_den = where(abs(del_off)-sigcut*abs(sig_off) gt 0)


;The difference quoted as a n-sigma uncertainty
;sig_den = abs(del_den[chk_den])/den_unc[chk_den]

;Use the Chauvenet's Criterion where sigma and <mu> come from the mission long comparision with WIND
;;;med_off = 6.3/100.
;;;;sig_off = 26.8/100.
;;;bad_tst = n_elements(del_den)*erfc(abs(del_den-med_off)/sig_off)
;;;bad_den = where(bad_tst lt 0.500000) ; Assume a gaussian distribution


print,n_elements(bad_den)

;where the measured uncertainty is 5 time greater than the DSCOVR WIND difference
;bad_den = where(sig_den gt sigcut)

;create array groups for areas where the density is 5 times greater than the unc.
;for an extended period (npix)
if n_elements(size(bad_den)) gt 3 then windcut = grp_obs(bad_den,npix) else windcut = -9999.0
if n_elements(size(windcut)) gt 3 then windcut = chk_den[temporary(windcut)] else windcut = -9999.0

;compute the array values to return to the main program
;;if n_elements(size(ind_grp)) gt 3 then windcut = bad_den[ind_grp] else windcut = -9999.0
;;
;;;print information to make sure its working properly (hint it is not)
;if n_elements(size(windcut)) gt 3 then print,sig_den[windcut],npix,sigcut

;check to see the median value inside the block of values

return,windcut
end

;--------------------------------------------------
;
;USAGE
;replace = sig_flag(y,yerr,npix=npix,sigcut=sigcut,tol=tol)
;
;COMMENTS
;Sends flag for values which signficantly differ from the running median
;--------------------------------------------------
function sig_flag,y,yerr,medar=medar,dmeda=dmeda,npix=npix,sigcut=sigcut
if keyword_set(npix) then npix = npix else npix = 1
if keyword_set(sigcut) then sigcut=sigcut else sigcut = 20.00;0.75

running_med,y,medar,dmeda,npix=npix


;find points sigcut away from the median
offset = abs(y-medar)/yerr
higsig = where((offset gt sigcut) and (abs(y-medar) gt 100),cnt)

return,higsig
end


;--------------------------------------------------
;
;This SUBROUTINE defines a idl strcture template to keep track of idl versioning
;
;USAGE
;file_version_temp
;
;--------------------------------------------------
pro file_version_temp

tmp = {file_data, $
       obs_year:fix(0),$
       obs_doy:fix(0),$
       idl_dir:'',$
       cdf_dir:'',$
       idl_file:'',$
       cdf_file:'',$
       idl_creation_date:'',$
       cdf_creation_date:'',$
       idl_creation_date_jd:double(0.0),$
       cdf_creation_date_jd:double(0.0)$
      }


end

;--------------------------------------------------
; SUBROUTINE for coverting fraction of a day (fd) into hour (hh), minute (mm), and second (ss) variables
;USAGE
;frac_hhmmss,fd,hh,mm,ss
;--------------------------------------------------
pro frac_hhmmss,fd,hh,mm,ss


fh = 24.*fd; day fraction in hours
hh = fix(fh) ; number of hours into day

fm = 60.*(fh-hh) ; remaining fraction in minutes
mm = fix(fm)     ; number of minutes into hour

fs = 60.*(fm-mm) ;remaining fraction in seconds
ss = fs     ;number of seconds into minute

end

;Function converts idl fill values to cdf fill values
function fill_bad,val,idlfil,cdffil
bad = where(val lt idlfil+1.)
if n_elements(size(bad)) gt 3 then val[bad]= cdffil
return,val
end

;--------------------------------------------------
;Function takes root idl structure and observed epoch and send fill value to bad values
;USAGE
;root = remove(remove_bad_times,root,obsepoch)
;--------------------------------------------------
function remove_bad_times,root,obsepoch

badfile = '../rejected_times/rejected_times.dat' ; file containing formatted rejected times


;create variables to hold year, doy, hour, min
syear = 0 & sdoy = 0 & sh = 0 & sm = 0
eyear = 0 & edoy = 0 & eh = 0 & em = 0

;create empty epoch arrays
strepoch = MAKE_ARRAY(1,/double)
endepoch = MAKE_ARRAY(1,/double)

;open file
openr,lun,badfile,/get_lun
;skip first two lines
skip_lun,lun,2,/LINES

;read format
rfmt = '(I4,1x,I03,1x,I02,1x,I02,1x,I4,1x,I03,1x,I02,1x,I02)'
;open bad file 
WHILE NOT EOF(lun) DO BEGIN & $
    readf, lun,format=rfmt,syear,sdoy,sh,sm,eyear,edoy,eh,em & $
    cdf_tt2000,tstrepoch,syear,0,sdoy,sh,sm,/compute_epoch & $;compute cdf epoch time for bad start
    cdf_tt2000,tendepoch,eyear,0,edoy,eh,em,/compute_epoch & $;compute cdf epoch time for bad end
    strepoch = [strepoch,tstrepoch] & $
    endepoch = [endepoch,tendepoch] & $
ENDWHILE
;close file
CLOSE,lun 
Free_lun,lun

; lop off padding 0
strepoch = strepoch[1:*]
endepoch = endepoch[1:*]



;loop over bad times and replace with fill value of -9999.0
fillval = -9999.0
for i=0,n_elements(strepoch)-1 do begin

;values in time range
    bad = where((obsepoch ge strepoch[i]) and (obsepoch le endepoch[i]))

    if n_elements(size(bad)) gt 3 then begin
;fill all arrays with bad values
        root.VX.data[bad] = fillval
        root.VY.data[bad] = fillval
        root.VZ.data[bad] = fillval
        
        root.VX.uncertainty[bad] = fillval
        root.VY.uncertainty[bad] = fillval
        root.VZ.uncertainty[bad] = fillval
    
        
        root.W.data[bad] = fillval 
        root.W.uncertainty[bad] = fillval 
        
        root.N.data[bad] = fillval
        root.N.uncertainty[bad] = fillval
    
        root.T.data[bad] = fillval
        root.T.uncertainty[bad] = fillval
    endif
endfor



return,root
end


;--------------------------------------------------
;ROUTINE for coverting idlsave file to CDF file
;Coverts idl save files to CDF format
;
;USAGE
;idl_dscovr_to_cdf,year,doy,version,archive=archive,skeleton=skeleton,filefmt=filefmt, $
;                  outfmt=outfmt,orchive=orchive,outdom=outdom
;--------------------------------------------------
pro idl_dscovr_to_cdf,year,doy,version,archive=archive,skeleton=skeleton,filefmt=filefmt,outfmt=outfmt,orchive=orchive,outdom=outdom

if keyword_set(archive) then archive=archive else archive='/crater/observatories/dscovr/plasmag/l2/idl/public_kp/1minute_corrected';get the default location of DSCOVR archive
archive = archive+'/'

if keyword_set(orchive) then orchive=orchive else orchive='/crater/observatories/dscovr/plasmag/l2/cdf/public_kp/1minute_corrected';set output cdf location
orchive = orchive+'/'

if keyword_set(skeleton) then skeleton=skeleton else skeleton='../skeleton/dscovr_h1_fc_0000000_v01.cdf';read the default cdf skeleton if none supplied

if keyword_set(filefmt) then filefmt=filefmt else filefmt='("dsc_fc_advkp_1minute_corrected_",I4,"_",I03,".idl")' ;default file format is for 1 minute cadence data

;if keyword_set(outfmt) then outfmt=outfmt else outfmt = '("dscovr_h1_plasmag_v01_",I4,I03,I02,I02,I02,".cdf")';old format when 1 file per observation
if keyword_set(outfmt) then outfmt=outfmt else outfmt = '("dscovr_h1_fc_",I4,I02,I02,"_v",I02,".cdf")' ;updated version by save file Prchlik J. (2017/03/16)


;base file name
bfil =string([year,doy],format=filefmt);get file name based on doy and year
file = archive+bfil



restore,file ;restore IDL save file

oroot = DFC_KP_1MIN
dayfrac = oroot.time.data-doy ;covert day time data into fractions of day




;loop over all day fractions and write to cdf file
;for i=0, 2 do begin;n_elements(dayfrac)-1 do begin
fd = dayfrac
frac_hhmmss,fd,hh,mm,ss ;get HH,MM,SS from fraction of day
; out_cdf = string([year,doy,hh,mm,round(ss),version],format=outfmt) ;create filename from day fraction

adoy = doy+intarr(n_elements(fd))
ayear = year+intarr(n_elements(fd))




cdf_tt2000,obsepoch,ayear,intarr(n_elements(fd)),adoy,hh,mm,round(ss),/compute_epoch ;compute cdf epoch time



;day fraction in milliseconds for time_pb5
dayms = round(fd*24.*3600.*1000.);day*(hour/day)*(s/hour)*(ms/s)


;fill value if times are determined to be bad
root = remove_bad_times(oroot,obsepoch) ; program to find and remove bad times

;OUTPUT format in YYYYMMDD format
cdf_tt2000,obsepoch,oyear,omon,odom,ohh,omm,oss,/break;create month and dom variables
epochfmt = '(I4,"-",I02,"-",I02,"T",I02,":",I02,":",I02)'
mepoch  = string([year,omon,odom,ohh,omm,oss],format=epochfmt)



;Create file format
out_cdf = string([year,fix(omon[0]),fix(odom[0]),version],format=outfmt) ;create filename from day fraction file per day format    
;read cdf skeleton
buf1 = read_master_cdf(skeleton,orchive+out_cdf) ;create cdf file from master skeleton

;time in PB5 format
pb5t = transpose([[year+intarr(n_elements(dayms))],[doy+intarr(n_elements(dayms))],[dayms]])


vgse = transpose([[root.VX.data],[root.VY.data],[root.VZ.data]]) ; Velocity in GSE 
ugse = transpose([[root.VX.uncertainty],[root.VY.uncertainty],[root.VZ.uncertainty]]) ;Uncertainty in velocity 

;most probable speed
mps = root.W.data ; most probable speed
dmps= root.W.uncertainty ; most probable speed


;Proton density
inp = root.N.data
idnp= root.N.uncertainty

;Wind temperature
temp = root.T.data
utemp= root.T.uncertainty




;-----------------------------------------------------
;bit/byte table (2017/06/23 J. Prchlik)
;     s | 
; f   p | b  n
; i d i | y  o
; l e k | t  t
; l n e | e  e
; 0 0 0 | 0  good
; 0 0 1 | 1  spiked in DSCOVR
; 0 1 0 | 2  prolongged sigma density mismatch WIND/DSCOVR 
; 1 0 0 | 3  data drop out (filled)
;-------| 4-7 Cannot duplicate flags
;
;
;-----------------------------------------------------
;data quailty flag
dqf_val = intarr(n_elements(temp));good (0) by default

;set up variables
num_rec = n_elements(obsepoch)
if keyword_set(outdom) then epoch = strarr(num_rec) else epoch = double(num_rec) ; Prchlik J. add keyword to output Epoch in dom format
epoch_old = double(num_rec)
time_pb5 = L64INDGEN(3,num_rec)
dqf = uint(num_rec)
v_gse = fltarr(3,num_rec)
v_gse_delta = fltarr(3,num_rec)
thermal_spd = float(num_rec)
thermal_spd_delta = float(num_rec)
np = float(num_rec)
np_delta = float(num_rec)
thermal_temp = float(num_rec)
thermal_temp_delta = float(num_rec)
;replace bad values with value in fillval for CDF
floatfil = -1.0e31
givenfil = -9999.0

;Replace IDL fill with CDF fill for velocity
for i=0,2 do begin
    vgse[i,*] = fill_bad(vgse[i,*],givenfil,floatfil)
    ugse[i,*] = fill_bad(ugse[i,*],givenfil,floatfil)
endfor


;fill most probable speed
mps = fill_bad(mps,givenfil,floatfil)
dmps = fill_bad(dmps,givenfil,floatfil)

;density fill values
inp = fill_bad(inp,givenfil,floatfil)
idnp = fill_bad(idnp,givenfil,floatfil)

;fill thermal temperature
temp = fill_bad(temp,givenfil,floatfil)
utemp = fill_bad(utemp,givenfil,floatfil)

;raise bad data flag for all points with 
badv = where((vgse[0,*] lt -9998.) or (vgse[1,*] lt -9998.) or $
             (vgse[2,*] lt -9998.) or $
             (ugse[0,*] lt -9998.) or (ugse[1,*] lt -9998.) or $
             (ugse[2,*] lt -9998.) or $
             (mps lt -9998.) or (inp lt -9998.) or (temp lt -9998.) or $
             (dmps lt -9998.) or (idnp lt -9998.) or (utemp lt -9998.) $
           )

;set bad dqf_val where bad points exist
if n_elements(size(badv)) gt 3 then dqf_val[badv] = 3

;check for Vx values 5 sigma away from the median
user_check = sig_flag(root.VX.data,root.VX.uncertainty,sigcut=5,npix=2)
if n_elements(size(user_check)) gt 3 then dqf_val[user_check] = 1

;check for density values 5 sigma away from WIND for more than 25 minutes
dens_check = den_flag(root.N.data,root.N.uncertainty,root.time.data,year,sigcut=3,npix=25)
if n_elements(size(dens_check)) gt 3 then dqf_val[dens_check] = 2

;fix for solar wind's aberration in Y component
;solab = 29.78 ;km/s
vgse[1,*] = vgse[1,*];-solab; Moved to advanced dsc_advanced_kp_fit_v2 (2017/04/07 J. Prchlik)


; put values into variables
epoch_old = obsepoch
if keyword_set(outdom) then epoch=mepoch else epoch = obsepoch ; Prchlik J. 2017/03/09 out put in YYYY/MM/DD format
time_pb5(0,*) = pb5t[0,*]
time_pb5(1,*) = pb5t[1,*]
time_pb5(2,*) = pb5t[2,*]
dqf = dqf_val
v_gse(0,*) = vgse[0,*]
v_gse(1,*) = vgse[1,*]
v_gse(2,*) = vgse[2,*]
v_gse_delta(0,*) = ugse[0,*]
v_gse_delta(1,*) = ugse[1,*]
v_gse_delta(2,*) = ugse[2,*]
thermal_spd =  mps
thermal_spd_delta = dmps
np = inp
np_delta = idnp
thermal_temp = temp
thermal_temp_delta = utemp




;set all data buffers
*buf1.Epoch.data             = epoch                ;epoch data
*buf1.Time_PB5.data          = time_pb5             ;PB5 time format
;    *buf1.FORMAT_TIME.data       = 'test'
;    *buf1.unit_time.data         = 'test'
*buf1.DQF.data               = dqf                  ;data quality flag
*buf1.V_GSE.data             = v_gse                ;velocity in GSE
;*buf1.cartesian.data         = 'x,y,z'
;*buf1.label_V_GSE.data       = 'km/s'
*buf1.V_GSE_DELTA.data       = v_gse_delta          ;uncertainty in V_GSE
;*buf1.label_V_GSE_DELTA.data = 'km/s'
*buf1.THERMAL_SPD.data       = thermal_spd          ;most probable speed
*buf1.THERMAL_SPD_DELTA.data = thermal_spd_delta    ;uncertainty in most probable speed
*buf1.Np.data                = np                   ;proton density
*buf1.Np_DELTA.data          = np_delta             ;uncertainty in proton density
*buf1.THERMAL_TEMP.data      = thermal_temp         ;proton wind temperature
*buf1.THERMAL_TEMP_DELTA.data= thermal_temp_delta   ;uncertainty in wind temperature



res1 = write_data_to_cdf(orchive+out_cdf,buf1,/debug) ;write output data to file

red_cdf = strmid(out_cdf,0,strlen(out_cdf)-4)
test = read_mycdf("Np,Np_DELTA",orchive+red_cdf,/all,/NODATASTRUCT) ; test to see if cdf write worked

print,orchive+red_cdf


;Add versioning information to idlsave file
savefmt = '("dscovr_file_status_v",I02,".idl")'
savefil = string([version],format=savefmt)
version_save = file_test(savefil)


;Get string of date and time of creation for idl and cdf files
SPAWN,"date -r "+archive+bfil+' +"%Y/%m/%d %H:%M:%S"',idl_date_str 
SPAWN,"date -r "+orchive+out_cdf+' +"%Y/%m/%d %H:%M:%S"',cdf_date_str 

;Get julian day for each string time
idl_date_jul = double(JULDAY(fix(strmid(idl_date_str,5,2)),fix(strmid(idl_date_str,8,2)),fix(strmid(idl_date_str,0,4)),fix(strmid(idl_date_str,11,2)),fix(strmid(idl_date_str,14,2)),fix(strmid(idl_date_str,17,2))))
cdf_date_jul = double(JULDAY(fix(strmid(cdf_date_str,5,2)),fix(strmid(cdf_date_str,8,2)),fix(strmid(cdf_date_str,0,4)),fix(strmid(cdf_date_str,11,2)),fix(strmid(cdf_date_str,14,2)),fix(strmid(cdf_date_str,17,2))))



tmp_s= {file_data, $
       obs_year:year,$
       obs_doy:doy,$
       idl_dir:archive,$
       cdf_dir:orchive,$
       idl_file:bfil,$
       cdf_file:out_cdf,$
       idl_creation_date:idl_date_str,$
       cdf_creation_date:cdf_date_str,$
       idl_creation_date_jd:idl_date_jul,$
       cdf_creation_date_jd:idl_date_jul $
      }


;if version file already exists concat structure else create new
if version_save eq 1 then begin
    restore,savefil,/RELAXED_STRUCTURE_ASSIGNMENT 
    v_struct = [v_struct,tmp_s]
    
endif else v_struct = tmp_s

;save the updated savefile
save,v_struct,filename=savefil

          

end 



