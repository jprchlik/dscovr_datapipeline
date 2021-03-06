
;--------------------------------------------------
;Uses 
;    load_plasma.pro
;Need to run the following if not running from dsc_advanced_kp_fit_v2.pro
;     @'/crater/utilities/idl/mike/idlstartup' 
;
;
;USAGE
;compare_wind_dscovr,year,doy,span=span,filefmt=filefmt,archive=archive
;
;COMMENTS
;    Compares WIND and DSCOVR data for a given day and plot puts plots with 
;    WIND and DSCOVR overplotted (_ontop_)
;    and DSCOVR-WIND (_delta_) in ../out_plots/
;--------------------------------------------------
;

;--------------------------------------------------
;USAGE
;    error_bars,xval,yval,err,p_err_x,p_err_y
;
;COMMENTS
;    creates error bars for plots 
;    Allow the syntax claims only Y errors you can flop x and y to plot x errors
;--------------------------------------------------
pro error_bars,xval,yval,err,p_err_x,p_err_y

p_err_x = [xval,xval]
p_err_y = [abs(err),-abs(err)]+yval


end

;===================================================
;
;Return newest version of file
;
;USAGE
;new_file = new_ver(archive,file)
;
;COMMENTS
;   Finds the most recent version of file to return
;===================================================
function new_ver,archive,file


new_file = file_search(archive+file)
if n_elements(size(new_file)) gt 3 then new_file = new_file[n_elements(new_file)-1]


return,new_file
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


for i=npix,sizer-npix-1 do begin
    medvl = y[i-npix:i+npix]
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
;USAGE
;dtime = chi_min_time(wdoy,wd,wx,wy,wz,ddoy,dd,dx,dy,dz,ddd,ddx,ddy,ddz,span=span,samp=samp
;
;COMMENTS
;    Function which spans +/- 10 minutes to find Chi^2 min
;    wdoy  = WIND day of year
;    wd = WIND-density 
;    wx = WIND-Vx      
;    wy = WIND-Vy      
;    wz = WIND-Vz      
;    ddoy  = DSCOVR day of year
;    dd = DSCOVR-density
;    dx = DSCOVR-Vx
;    dy = DSCOVR-Vy
;    dz = DSCOVR-Vz
;    ddd = DSCOVR uncertainty in density
;    ddx = DSCOVR uncertainty in Vx
;    ddy = DSCOVR uncertainty in Vy
;    ddz = DSCOVR uncertainty in Vz
;    span = time in minutes to minimize
;    samp = sampling fequence in minutes
;--------------------------------------------------
function chi_min_time_2,wdoy,wd,wx,wy,wz,ddoy,dd,dx,dy,dz,ddd,ddx,ddy,ddz,span=span,samp=samp
if keyword_set(span) then span = span else span = 10 ;span to loop +/- in minutes
if keyword_set(samp) then samp = samp else samp = 0.50 ; sampling fequency in minutes

span = span/60./24. ;turn into fraction of a day
samp = samp/60./24. ;turn into fraction of a day

;array of offset to loop over
nbins = round(2.*span/samp)+1
tvals = findgen(nbins)*samp-span
cvals = fltarr(nbins) ; array of zeros to store chi^2 vals


;loop over points for chi2 sampling
for i=0,nbins-1 do begin
    ;create spline from wind data assuming small sigma
    sdd = spline(wdoy,wd,ddoy+tvals[i])-dd
    sdx = spline(wdoy,wx,ddoy+tvals[i])-dx
    sdy = spline(wdoy,wy,ddoy+tvals[i])-dy
    sdz = spline(wdoy,wz,ddoy+tvals[i])-dz
    chi2 = total(sdd^2/ddd^2)+total(sdx^2/ddx^2)+total(sdy^2/ddy^2)+total(sdz^2/ddx^2)
;store chi2 val
    cvals[i] = chi2
endfor

;fine sample chi2 distirbution to find minimum value
fsamp = samp/100.
;array of offset to loop over
fnbins = round(2.*span/fsamp)+1
ftvals = findgen(fnbins)*fsamp-span
fcvals = spline(tvals,cvals,ftvals) ; array of of chi^2 vals

;minimum chi^2 
chim = where(fcvals eq min(fcvals))
if n_elements(size(chim)) gt 3 then begin ;move the time as little as possible
    dtime = ftvals[chim]
    minmo = where(abs(dtime) eq min(abs(dtime))) ; move the minimum amount
    dtime = dtime[minmo]
endif else begin 
    dtime = 0. ; return 0 if no min??
endelse
;print,'Chi^2 min value',strcompress(min(fcvals),/remove_all)
;print,'Chi^2 time offset',strcompress(dtime*24.*60.),'min'


return,dtime
end

;--------------------------------------------------
;MAIN ROUTINE
;USAGE
;compare_wind_dscovr,year,doy,version,span=span,filefmt=filefmt,archive=archive
;
;COMMENTS
;    Compares WIND and DSCOVR data for a given day and plot puts plots with 
;    WIND and DSCOVR overplotted (_ontop_)
;    and DSCOVR-WIND (_delta_) in ../out_plots/
;    year and day of year (doy) are integers 
;    span in days to compare WIND and DSCOVR (default = 1)
;    filefmt is the format of the idl save file you need to restore (default = '("dscovr_h1_fc_",I4,I02,I02,"_v",I02)')
;    archive is the location of the cdf files archive (default = '/crater/observatories/dscovr/plasmag/l2/cdf/public_kp/1minute_corrected)
;    ffmt is the file format of the output png file with plasma dscovr plotted ontop of wind
;    dfmt is the file format of the output png file with dscovr wind plasma difference plotted 
;    mfmt is the file format of the output png file with magnetic field dscovr plotted ontop of wind
;
;--------------------------------------------------
pro compare_wind_dscovr,year,doy,version,span=span,filefmt=filefmt,archive=archive,ffmt=ffmt,mfmt=mfmt,dfmt=dfmt

if keyword_set(archive) then archive=archive else archive='/crater/observatories/dscovr/plasmag/l2/cdf/public_kp/1minute_corrected';get the default location of DSCOVR archive
archive= archive+'/'
if keyword_set(filefmt) then filefmt=filefmt else filefmt= '("dscovr_h1_fc_",I4,I02,I02,"_v",I02)';default file format is for 1 minute cadence data
if keyword_set(span) then span=fix(span) else span = 1
if keyword_set(ffmt) then ffmt = ffmt else ffmt = '("../out_plots/compare_wind_dsc_ontop_",I4,"_",I03,".png")'
if keyword_set(dfmt) then dfmt = dfmt else dfmt = '("../out_plots/compare_wind_dsc_delta_",I4,"_",I03,".png")'
if keyword_set(mfmt) then mfmt = mfmt else mfmt = '("../out_plots/compare_wind_dsc_ontop_mag_",I4,"_",I03,".png")'

plotsym,0,/fill

;Convert DOY to YYYYMMDD format
byear = JULDAY(01,01,year) ;First day of year
bdoyj = byear+doy-1
CALDAT,bdoyj,mon,dom,yyyy

;file = archive+string([year,doy],format=filefmt);get file name based on doy and year
file = archive+string([yyyy,mon,dom,version],format=filefmt)

;WIND data
load_plasma,year,doy,doy+span,/wind,vx=wx,vy=wy,vz=wz,t=wt,den=wd,doy=wdoy,ns=wew
;ACE data from 2008 (i.e. a leap year)
lp = year
load_plasma,lp,doy,doy+span,/ace,vx=ax,vy=ay,vz=az,t=at,den=ad,doy=adoy,ns=aew

;convert thermal speed to mp speed
;constants
amu = double(1.660538921e-27);amu
mp  = 1.00727647*amu
kb  = double(1.3906488e-23);boltzmann's constant

;convert temperature to most probable speed
wew = wt
aew = at

;restore,file ;restore idl save file
root = read_mycdf("",file,/all)
;root = DFC_KP_1MIN

;ddoy = root.time.data ; Day values
epoch = root.epoch.dat
cdf_epoch,epoch,yr,mo,dm,hr,mn,sc,milli,/break

;VX
dx = root.v_gse.dat[0,*]
ddx= root.v_gse_delta.dat[0,*]


;VY
dy = root.v_gse.dat[1,*]
ddy= root.v_gse_delta.dat[1,*]


;VZ
dz = root.v_gse.dat[2,*]
ddz= root.v_gse_delta.dat[2,*]

;Proton density
dd = root.NP.dat
ddd= root.NP_DELTA.dat

;Proton thermal speed 
dew = 1.E-3*sqrt(2.*kb/mp*root.THERMAL_TEMP.dat)
ddew= 1.E-3*sqrt(2.*kb/mp/root.THERMAL_TEMP.dat)*root.THERMAL_TEMP_DELTA.dat*0.5
;ddew=0.5*(root.THERMAL_TEMP_DELTA.dat/root.THERMAL_TEMP.dat) errbar sanity check

;data quality flag
dqf = root.dqf.dat




;============================================================
; Find Mag. data if avialable Section 1.1,1
;============================= ===============================
;Set up archive and file information for trajectory and magnetic field archives
if keyword_set(marchive) then marchive=marchive else marchive='/crater/observatories/dscovr/plasmag/mag';set output cdf location
marchive = marchive+'/'
if keyword_set(mfilefmt) then mfilefmt=mfilefmt else mfilefmt = '(I4,"/dscovr_h0_mag_",I4,I02,I02,"_v*.cdf")'

;find path for cdf file
;mcdf = archive+string([year,year,mon,dom,v1],format=mfilefmt)

;find path for orbital file
mcdf = string([year,year,mon,dom],format=mfilefmt)
mag_cdf = new_ver(marchive,mcdf)

;make sure files exsits
test_mag = file_test(mag_cdf)

;if file exits read it n
if (test_mag eq 1) then mag_cdf = read_mycdf("",mag_cdf,/all)


;============================================================
; Set up plots Section 1.2
;============================================================

;device,decompsed=0
;device,retain=2 ;prevent from plotting empty windows

;image size

set_plot,'Z'
device,set_resolution=[1500,1000],decomposed=0,set_pixel_depth=24
loadct,12
;device,xsize=8.5,ysize=11.0,/inches
!p.thick=4
!x.thick=3
!y.thick=3
;!y.margin=[0,0]
;!x.margin=[15,0]
;!P.MULTI = [0,1,4]

ecol = 224 ;errorbar color
wcol = 200 ; red wind color
dcol = 0   ;black dscovr color
acol = 100 ; ACE color


;covert doy to JULDATES
;jddoy = double(julday(1,1,year,0,0,0)+(ddoy-1)) ; DSCOVR
jddoy = double(julday(mo,dm,yr,hr,mn,sc))
ddoy = jddoy-julday(1,1,year,0,0,0)+1.; covert to doy
jwdoy = double(julday(1,1,year,0,0,0)+(wdoy-1)) ; WIND  

;setup format for plots
dummy = label_date(date_format=['%H:%I'])
;set xrange for plots
twoh = double(0.083)
xrange = [min(jddoy)-twoh,max(jddoy)+0.084333] ;pad 2hours
!X.range = xrange
!X.style = 1
!X.ticks = 14
!X.minor = 4
!X.ticklen = 0.05
!X.ticks = 14
!Y.ticks = 4
xticks = ['22:00','00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','00:00','02:00']


xlabfmt = '("DOY=",I03,"/",I4)'


;set up arrays for spanning each plot range
yb0 = .08
yt5 = .95
ysp = yt5-yb0 ;obsevation span
psz = ysp/5.  ; size of each oplot


;set up plot positions for each plot
plot1 = [.1,yb0+4.*psz,.95,yb0+5.*psz]
plot2 = [.1,yb0+3.*psz,.95,yb0+4.*psz]
plot3 = [.1,yb0+2.*psz,.95,yb0+3.*psz]
plot4 = [.1,yb0+1.*psz,.95,yb0+2.*psz]
plot5 = [.1,yb0+0.*psz,.95,yb0+1.*psz]


;============================================================
;PLOT SECTION 1.2.1 
;Overplot WIND AND DSCVOR
;============================================================


;Plot proton density
;where the density if off by n sigma
dn = where(dqf eq 2,dncnt)
plot,jwdoy,wd,psym=6,color=0,ytitle='Density [cm^-3]',/nodata,background=255,charsize=2,$
     font=1,charthick=3,position=plot1,xtickformat="(A1)"
for i=0,n_elements(jddoy)-1 do begin
    error_bars,jddoy[i],dd[i],ddd[i],ex,ey
    if dd[i] gt -9998. then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
endfor
oplot,jwdoy,wd,psym=6,color=wcol
oplot,jddoy,dd,psym=8,color=dcol
if dncnt gt 0 then oplot,jddoy[dn],dd[dn],psym=8,color=120,s=2.0
;oplot,adoy,ad,psym=5,color=acol

;Plot proton thermal speed
plot,jwdoy,wew,psym=6,color=0,ytitle='Th. Speed [km/s]',/nodata,background=255,charsize=2,$
     font=1,charthick=3,position=plot2,/NOERASE,xtickformat="(A1)"
for i=0,n_elements(jddoy)-1 do begin
    error_bars,jddoy[i],dew[i],ddew[i],ex,ey
    if dew[i] gt -9998. then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
endfor
oplot,jwdoy,wew,psym=6,color=wcol
oplot,jddoy,dew,psym=8,color=dcol
if dncnt gt 0 then oplot,jddoy[dn],dew[dn],psym=8,color=120,s=2.0
;oplot,adoy,ad,psym=5,color=acol


;Plot Vx
ui = where(dqf eq 1,uicnt)
plot,jwdoy,wx,psym=6,color=0,ytitle='Vx [km/s]',/nodata,background=255,charsize=2,$
     font=1,charthick=3,position=plot3,/NOERASE,xtickformat="(A1)"
;plot runnint median
running_med,dx,medar,dmeda,npix=2
oplot,jddoy,medar,psym=10,color=120,linestyle=0,thick=.2
for i=0,n_elements(jddoy)-1 do begin
    error_bars,jddoy[i],dx[i],ddx[i],ex,ey
    if dx[i] gt -9998. then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
endfor
oplot,jwdoy,wx,psym=6,color=wcol
oplot,jddoy,dx,psym=8,color=dcol
if uicnt gt 0 then oplot,jddoy[ui],dx[ui],psym=8,color=100,s=2.0
;oplot,adoy,ax,psym=5,color=acol

;Plot Vy
plot,jwdoy,wy,psym=6,color=0,ytitle='Vy [km/s]',/nodata,background=255,charsize=2,$
     font=1,charthick=3,position=plot4,/NOERASE,xtickformat="(A1)"
for i=0,n_elements(jddoy)-1 do begin
    error_bars,jddoy[i],dy[i],ddy[i],ex,ey
    if dy[i] gt -9998. then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
endfor
oplot,jwdoy,wy,psym=6,color=wcol
oplot,jddoy,dy,psym=8,color=dcol
;oplot,adoy,ay,psym=5,color=acol

;Plot Vz
plot,jwdoy,wz,psym=6,color=0,ytitle='Vz [km/s]',xtitle=string([doy,year],format=xlabfmt),/nodata,background=255,charsize=2,$
     font=1,charthick=3,position=plot5,/NOERASE,xtickformat='label_date'
for i=0,n_elements(jddoy)-1 do begin
    error_bars,jddoy[i],dz[i],ddz[i],ex,ey
    if dz[i] gt -9998. then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
endfor
oplot,jwdoy,wz,psym=6,color=wcol
oplot,jddoy,dz,psym=8,color=dcol
;oplot,adoy,az,psym=5,color=acol


fname = string([year,doy],format=ffmt)
TVLCT,r,g,b,/Get
write_png,fname,tvrd(/true),r,g,b
;write_png,fname,tvrd(/true)
;device,/close


;============================================================
;Section 1.2.2
;Find DSCOVR magnetic field data
;============================================================

;If mag cdf file exists create plots
if (test_mag eq 1) then begin

    ;set up arrays for spanning each plot range
    yb0 = .08
    yt5 = .95
    ysp = yt5-yb0 ;obsevation span
    psz = ysp/4.  ; size of each oplot
    
    
    ;set up plot positions for each plot
    plot2a = [.1,yb0+3.*psz,.95,yb0+4.*psz]
    plot3a = [.1,yb0+2.*psz,.95,yb0+3.*psz]
    plot4a = [.1,yb0+1.*psz,.95,yb0+2.*psz]
    plot5a = [.1,yb0+0.*psz,.95,yb0+1.*psz]

    ;Good quailty magnometer flag
    good_mag = where(mag_cdf.flag1.dat eq 0)

    ;convert cdf epoch to julian days
    mag_time = CDF_EPOCH_TOJULDAYS(mag_cdf.epoch1.dat)

    plot,mag_time[good_mag],mag_cdf.b1f1.dat[good_mag],psym=6,color=0, $
     ytitle='|B| [nT]',background=255,charsize=2,$
     font=1,charthick=3,position=plot2a,xtickformat="(A1)"
    ;for i=0L,n_elements(good_mag)-1 do begin
    ;    error_bars,mag_time[good_mag[i]],mag_cdf.b1f1.dat[good_mag[i]],mag_cdf.b1sdf1.dat[good_mag[i]],ex,ey
    ;    oplot,ex,ey,color=ecol,linestyle=0,thick=.2
    ;endfor

    plot,mag_time[good_mag],mag_cdf.b1gse.dat[0,good_mag],psym=6,color=0, $
     ytitle='Bx [nT]',background=255,charsize=2,$
     font=1,charthick=3,position=plot3a,/NOERASE,xtickformat="(A1)"
    ;for i=0L,n_elements(good_mag)-1 do begin
    ;    j = good_mag[i]
    ;    error_bars,mag_time[j],mag_cdf.b1gse.dat[0,j],mag_cdf.b1sdgse.dat[0,j],ex,ey
    ;    oplot,ex,ey,color=ecol,linestyle=0,thick=.2
    ;endfor

    plot,mag_time[good_mag],mag_cdf.b1gse.dat[1,good_mag],psym=6,color=0, $
     ytitle='By [nT]',background=255,charsize=2,$
     font=1,charthick=3,position=plot4a,/NOERASE,xtickformat="(A1)"
    ;for i=0L,n_elements(good_mag)-1 do begin
    ;    j = good_mag[i]
    ;    error_bars,mag_time[j],mag_cdf.b1gse.dat[1,j],mag_cdf.b1sdgse.dat[1,j],ex,ey
    ;    oplot,ex,ey,color=ecol,linestyle=0,thick=.2
    ;endfor

    plot,mag_time[good_mag],mag_cdf.b1gse.dat[2,good_mag],psym=6,color=0, $
     ytitle='Bz [nT]',xtitle=string([doy,year],format=xlabfmt),background=255,charsize=2,$
     font=1,charthick=3,position=plot5a,/NOERASE,xtickformat='label_date'
    ;oplot,mag_time[good_mag],mag_cdf.b1gse.dat[2,good_mag],psym=6
    ;for i=0L,n_elements(good_mag)-1 do begin
    ;    j = good_mag[i]
    ;    error_bars,mag_time[j],mag_cdf.b1gse.dat[2,j],mag_cdf.b1sdgse.dat[2,j],ex,ey
    ;    oplot,ex,ey,color=ecol,linestyle=0,thick=.2
    ;endfor

    fname = string([year,doy],format=mfmt)
    TVLCT,r,g,b,/Get
    write_png,fname,tvrd(/true),r,g,b

endif
;============================================================
;Section 1.2.3
;Find DSCOVR-WIND offsets
;============================================================
;WIND data for +/- 1 day from DSCOVR span (allows movement data for correlation
load_plasma,year,doy-1,doy+span+1,/wind,vx=bwx,vy=bwy,vz=bwz,t=bwt,den=bwd,doy=bwdoy,ns=bwew
load_plasma,lp,doy,doy+span,/wind,vx=box,vy=boy,vz=boz,t=bot,den=bod,doy=bodoy,ns=boew ; Wind data from the same time as ACE
;convert temperature to most probable speed
bwew = bwt
boew = bot


;Find all no fill DSCOVR data points
nofill = where((dd gt -9998.) and (dx gt -9998.) and (dy gt -9998.) and (dz gt -9998.))
ferrs = fltarr(n_elements(nofill))+1
;minimize differences using time and errors
dtime = chi_min_time_2(bwdoy,bwd,bwx,bwy,bwz,ddoy[nofill],dd[nofill],dx[nofill],dy[nofill],dz[nofill],ddd[nofill],ddx[nofill],ddy[nofill],ddz[nofill],span=12);output in fraction of day of hear
dtime = float(dtime[0])

;Compare deltas
;DSCOVR
;create spline from wind data assuming small sigma
sdd = spline(bwdoy,bwd,ddoy+dtime)
sdx = spline(bwdoy,bwx,ddoy+dtime)
sdy = spline(bwdoy,bwy,ddoy+dtime)
sdz = spline(bwdoy,bwz,ddoy+dtime)
sdew= spline(bwdoy,bwew,ddoy+dtime)

;differences in interpolated wind and dscovr value
dwdd = dd-sdd
dwdx = dx-sdx
dwdy = dy-sdy
dwdz = dz-sdz
dwdw = dew-sdew

;ACE
;Find all no fill ACE data points
nofill = where((ad gt -998.) and (ax gt -998.) and (ay gt -998.) and (az gt -998.))
ace=1 ;overplot ACE
if n_elements(nofill) lt 4 then ace=0 
;minimize differences using time and errors
if ace eq 1 then begin
    ferrs = fltarr(n_elements(nofill))+1
    dtime1 = chi_min_time_2(bodoy,bod,box,boy,boz,adoy[nofill],ad[nofill],ax[nofill],ay[nofill],az[nofill],ferrs,ferrs,ferrs,ferrs);output in fraction of day of hear
    dtime1 = float(dtime1[0])
    ;create spline from wind data assuming small sigma
    sad = spline(bodoy,bod,adoy+dtime1)
    sax = spline(bodoy,box,adoy+dtime1)
    say = spline(bodoy,boy,adoy+dtime1)
    saz = spline(bodoy,boz,adoy+dtime1)
    saew= spline(bodoy,boew,ddoy+dtime)
    
    ;differences in interpolated wind and dscovr value
    dwad = ad-sad
    dwax = ax-sax
    dway = ay-say
    dwaz = az-saz
    dwaw = aew-saew
endif
    



;Find DSCOVR wind outliers
sigout = 5.

;covert doy to JULDATES
ddoy = double(julday(1,1,year,0,0,0)+(ddoy-1)) ; DSCOVR
wdoy = double(julday(1,1,year,0,0,0)+(wdoy-1)) ; WIND  
;setup add julday if ace is in date range
if ace eq 1 then adoy = double(julday(1,1,year)+(adoy-1))


;Plot proton density
use = where(dd gt -9998.,c_den)
if c_den gt 0 then begin
    plot,ddoy[use],dwdd[use],psym=6,color=0,ytitle='Density [cm^-3]',$
         title='(DSCOVR-WIND) '+strcompress(year,/remove_all)+'; Time offset = '+strcompress(dtime*24.*3600.)+'s',$
        /nodata,background=255,charsize=2,font=1,charthick=3,position=plot1,xtickformat="(A1)"
    for i=0,n_elements(ddoy)-1 do begin
        error_bars,ddoy[i],dwdd[i],ddd[i],ex,ey
        if ddd[i] gt -9998.0 then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
    endfor
    oplot,ddoy,dwdd,psym=8,color=dcol
    if ace eq 1 then oplot,adoy,dwad,psym=5,color=acol
   ;plot outliers in red
   out = where(abs(dwdd) gt sigout*ddd)
   if n_elements(out) gt 2 then oplot,ddoy[out],dwdd[out],psym=8,color=200

endif



;Plot proton thermal speed
use = where(dew gt -9998.,n_dew)
if n_dew gt 0 then begin
    plot,ddoy[use],dwdw[use],psym=6,color=0,ytitle='Th. Speed [km/s]',/NOERASE, $
        /nodata,background=255,charsize=2,font=1,charthick=3,position=plot2,xtickformat="(A1)"
    for i=0,n_elements(ddoy)-1 do begin
        error_bars,ddoy[i],dwdw[i],ddew[i],ex,ey
        if dew[i] gt -9998.0 then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
    endfor
    oplot,ddoy,dwdw,psym=8,color=dcol
    if ace eq 1 then oplot,adoy,dwaw,psym=5,color=acol
    
    ;plot outliers in red
    out = where(abs(dwdd) gt sigout*ddew)
    if n_elements(out) gt 2 then oplot,ddoy[out],dwdw[out],psym=8,color=200
endif



;Plot Vx
use = where(dx gt -9998.,n_dx)
if n_dx gt 0 then begin
    plot,ddoy[use],dwdx[use],psym=6,color=0,ytitle='Vx [km/s]',/nodata,background=255,charsize=2,$
         font=1,charthick=3,position=plot3,/NOERASE,xtickformat="(A1)"
    for i=0,n_elements(ddoy)-1 do begin
        error_bars,ddoy[i],dwdx[i],ddx[i],ex,ey
        if ddx[i] gt -9998.0 then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
    endfor
    oplot,ddoy,dwdx,psym=8,color=dcol
    if uicnt gt 0 then oplot,ddoy[ui],dwdx[ui],psym=8,color=100,s=2.0
    if ace eq 1 then oplot,adoy,dwax,psym=5,color=acol
    
    ;plot outliers in red
    out = where(abs(dwdx) gt 5.*ddx)
    if n_elements(out) gt 2 then oplot,ddoy[out],dwdx[out],psym=8,color=200
endif

;Plot Vy
use = where(dy gt -9998.,n_dy)
if n_dy gt 0 then begin
    plot,ddoy[use],dwdy[use],psym=6,color=0,ytitle='Vy [km/s]',/nodata,background=255,charsize=2,$
         font=1,charthick=3,position=plot4,/NOERASE,xtickformat="(A1)"
    for i=0,n_elements(ddoy)-1 do begin
        error_bars,ddoy[i],dwdy[i],ddy[i],ex,ey
        if ddy[i] gt -9998.0 then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
    endfor
    oplot,ddoy,dwdy,psym=8,color=dcol
    if ace eq 1 then oplot,adoy,dway,psym=5,color=acol
    
    ;plot outliers in red
    out = where(abs(dwdy) gt 5.*ddy)
    if n_elements(out) gt 2 then oplot,ddoy[out],dwdy[out],psym=8,color=200
endif

;Plot Vz
use = where(dz gt -9998.,n_dz)
if n_dz gt 0 then begin
    plot,ddoy[use],dwdz[use],psym=6,color=0,ytitle='Vz [km/s]',xtitle=string([doy,year],format=xlabfmt),/nodata,background=255,charsize=2,$
         font=1,charthick=3,position=plot5,/NOERASE,xtickformat='label_date'
    for i=0,n_elements(ddoy)-1 do begin
        error_bars,ddoy[i],dwdz[i],ddz[i],ex,ey
        if ddz[i] gt -9998.0 then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
    endfor
    oplot,ddoy,dwdz,psym=8,color=dcol
    if ace eq 1 then oplot,adoy,dwaz,psym=5,color=acol
    
    ;plot outliers in red
    out = where(abs(dwdz) gt 5.*ddz)
    if n_elements(out) gt 2 then oplot,ddoy[out],dwdz[out],psym=8,color=200
endif

fname = string([year,doy],format=dfmt)

TVLCT,r,g,b,/Get
write_png,fname,tvrd(/true),r,g,b

end
