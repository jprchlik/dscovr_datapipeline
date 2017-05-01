;----------------------------------------------------------------------------------------
;
;Uses
;
;USAGE 
;compare_wind_dscovr_ace,syear,sdoy,edoy=edoy,aceoff=aceoff,eyear=eyear
;
;
;----------------------------------------------------------------------------------------


;--------------------------------------------------
;Uses
;
;
;
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
function chi_min_time,wdoy,wx,ddoy,dx,ddx=ddx,span=span,samp=samp
if keyword_set(span) then span = span else span = 10 ;span to loop +/- in minutes
if keyword_set(samp) then samp = samp else samp = 0.50 ; sampling fequency in minutes
if keyword_set(ddx) then ddx = ddx else ddx = fltarr(n_elements(dx))+1.

span = span/60./24. ;turn into fraction of a day
samp = samp/60./24. ;turn into fraction of a day

;array of offset to loop over
nbins = round(2.*span/samp)+1
tvals = findgen(nbins)*samp-span
cvals = fltarr(nbins) ; array of zeros to store chi^2 vals


;loop over points for chi2 sampling
for i=0,nbins-1 do begin
    ;create spline from wind data assuming small sigma
;    sdd = spline(wdoy,wd,ddoy+tvals[i])-dd
    sdx = spline(wdoy,wx,ddoy+tvals[i])-dx
;    sdy = spline(wdoy,wy,ddoy+tvals[i])-dy
;    sdz = spline(wdoy,wz,ddoy+tvals[i])-dz
    chi2 = total(sdx^2/ddx^2)
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
print,'Chi^2 min value',strcompress(min(fcvals),/remove_all)
print,'Chi^2 time offset',strcompress(dtime*24.*60.),'min'


return,dtime
end


;----------------------------------------------------------------------------------------
;
;
;USAGE
;function time_offset_cor,wdoy,wvx,doy,vx,dvx=dvx
;
;
;----------------------------------------------------------------------------------------

function time_offset_cor,wdoy,wvx,doy,vx,dvx=dvx

if keyword_set(dvx) then dvx = dvx else dvx = fltarr(n_elements(vx))+1.

iday = fix(doy)
days = uniq(iday)
days = iday[days]
new_doy = fltarr(n_elements(doy))

wcut = where(wvx gt -9990.0)
;loop over all days in year
for i =0,n_elements(days) - 1 do begin 
    cut = where((iday eq days[i]) and (vx gt -9990.0))
    deldoy = chi_min_time(wdoy[wcut],wvx[wcut],doy[cut],vx[cut],ddx=dvx)
    new_doy[cut] = doy[cut]+deldoy;add time offset to new array

endfor

return,new_doy
end


;----------------------------------------------------------------------------------------
;
;
;USAGE
;format_hist,hist,bins,bsize
;
;
;----------------------------------------------------------------------------------------
pro format_hist,hist,bins,bsize
hist = [0,hist,0]
hist = 100*hist/float(total(hist))
bins = [bins[0]-bsize,bins,bins[n_elements(bins)-1]+bsize]
end

pro compare_wind_dscovr_ace,syear,sdoy,edoy=edoy,aceoff=aceoff,eyear=eyear

if keyword_set(aceoff) then aceoff=aceoff else aceoff = 11. ;offset for ace data in years

if keyword_set(edoy) then edoy=edoy else edoy = sdoy+1

if keyword_set(eyear) then eyear=eyear else eyear = syear

;array of year values to loop through
ryear = syear+findgen(eyear-syear+1)
ayear = ryear-aceoff ;array of ace offset years
;number of years
nyear = n_elements(ryear)


;loop over all years and store variables
for i=0,nyear-1 do begin

    if i eq nyear-1 then tedoy = edoy else tedoy = 365

    if i eq 0 then begin
        tsdoy = sdoy
        load_plasma,ryear[i],tsdoy,tedoy,/DSC,doy=ddoy,bx=dbx,by=dby,bz=dbz,bmag=dbmag,vmag=dvmag,vy=dvy,vx=dvx,vz=dvz,t=dt,den=dden,x=dx,y=dy,z=dz
        load_plasma,ryear[i],tsdoy,tedoy,/WIND,doy=wdoy,bx=wbx,by=wby,bz=wbz,bmag=wbmag,vmag=wvmag,vy=wvy,vx=wvx,vz=wvz,t=wt,den=wden,x=wx,y=wy,z=wz

;correct for orbital time offsets DSCOVR/WIND
        ddoy = time_offset_cor(wdoy,wvx,ddoy,dvx)

;do the same for ace
        load_plasma,ayear[i],tsdoy,tedoy,/ACE,doy=adoy,bx=abx,by=aby,bz=abz,bmag=abmag,vmag=avmag,vy=avy,vx=avx,vz=avz,t=at,den=aden,x=ax,y=ay,z=az
        load_plasma,ayear[i],tsdoy,tedoy,/WIND,doy=idoy,bx=ibx,by=iby,bz=ibz,bmag=ibmag,vmag=ivmag,vy=ivy,vx=ivx,vz=ivz,t=it,den=iden,x=ix,y=iy,z=iz
;correct for orbital time offsets ACE/WIND
        adoy = time_offset_cor(wdoy,wvx,adoy,avx)

    endif else begin
         tsdoy = 1
         load_plasma,ryear[i],tsdoy,tedoy,/DSC,doy=tddoy,bx=tdbx,by=tdby,bz=tdbz,bmag=tdbmag,vmag=tdvmag,vy=tdvy,vx=tdvx,vz=tdvz,t=tdt,den=tdden,x=tdx,y=tdy,z=tdz
         load_plasma,ryear[i],tsdoy,tedoy,/WIND,doy=twdoy,bx=twbx,by=twby,bz=twbz,bmag=twbmag,vmag=twvmag,vy=twvy,vx=twvx,vz=twvz,t=twt,den=twden,x=twx,y=twy,z=twz
;correct for orbital time offsets DSCOVR/WIND
         tddoy = time_offset_cor(wdoy,wvx,tddoy,tdvx)

;do the same for ace
         load_plasma,ayear[i],tsdoy,tedoy,/ACE,doy=tadoy,bx=tabx,by=taby,bz=tabz,bmag=tabmag,vmag=tavmag,vy=tavy,vx=tavx,vz=tavz,t=tat,den=taden,x=tax,y=tay,z=taz
         load_plasma,ayear[i],tsdoy,tedoy,/WIND,doy=tidoy,bx=tibx,by=tiby,bz=tibz,bmag=tibmag,vmag=tivmag,vy=tivy,vx=tivx,vz=tivz,t=tit,den=tiden,x=tix,y=tiy,z=tiz
;correct for orbital time offsets ACE/WIND
         tadoy = time_offset_cor(wdoy,wvx,tadoy,tavx)



;Put dscovr data in large arrays
         ddoy = [ddoy ,tddoy+365.0*i]
         dbx  = [dbx  ,tdbx]
         dby  = [dby  ,tdby]
         dbz  = [dbz  ,tdbz]
         dbmag= [dbmag,tdbmag]
         dvmag= [dvmag,tdvmag]
         dvy  = [dvy  ,tdvy]
         dvx  = [dvx  ,tdvx]
         dvz  = [dvz  ,tdvz]
         dt   = [dt   ,tdt]
         dden = [dden ,tdden]
         dx   = [dx   ,tdx]
         dy   = [dy   ,tdy]
         dz   = [dz   ,tdz]

;Put wind data in large arrays
         wdoy = [wdoy ,twdoy+365.0*i]
         wbx  = [wbx  ,twbx]
         wby  = [wby  ,twby]
         wbz  = [wbz  ,twbz]
         wbmag= [wbmag,twbmag]
         wvmag= [wvmag,twvmag]
         wvy  = [wvy  ,twvy]
         wvx  = [wvx  ,twvx]
         wvz  = [wvz  ,twvz]
         wt   = [wt   ,twt]
         wden = [wden ,twden]
         wx   = [wx   ,twx]
         wy   = [wy   ,twy]
         wz   = [wz   ,twz]

;Put ace data in large arrays
         adoy = [adoy ,tadoy+365.0*i]
         abx  = [abx  ,tabx]
         aby  = [aby  ,taby]
         abz  = [abz  ,tabz]
         abmag= [abmag,tabmag]
         avmag= [avmag,tavmag]
         avy  = [avy  ,tavy]
         avx  = [avx  ,tavx]
         avz  = [avz  ,tavz]
         at   = [at   ,tat]
         aden = [aden ,taden]
         ax   = [ax   ,tax]
         ay   = [ay   ,tay]
         az   = [az   ,taz]

;Put ace-wind data in large arrays
         idoy = [idoy ,tidoy+365.0*i]
         ibx  = [ibx  ,tibx]
         iby  = [iby  ,tiby]
         ibz  = [ibz  ,tibz]
         ibmag= [ibmag,tibmag]
         ivmag= [ivmag,tivmag]
         ivy  = [ivy  ,tivy]
         ivx  = [ivx  ,tivx]
         ivz  = [ivz  ,tivz]
         it   = [it   ,tit]
         iden = [iden ,tiden]
         ix   = [ix   ,tix]
         iy   = [iy   ,tiy]
         iz   = [iz   ,tiz]




    endelse
endfor

fillval = -9990.0

;Remove bad information DSCOVR
dgood = where((dbx gt fillval) and (dby gt fillval) and (dbz gt fillval) and (dbmag gt fillval) and (dvmag gt fillval) $
    and (dvx gt fillval) and (dvy gt fillval) and (dvz gt fillval) and (dt gt fillval) and (dden gt fillval)) 
ddoy = ddoy [dgood]
dbx  = dbx  [dgood]
dby  = dby  [dgood]
dbz  = dbz  [dgood]
dbmag= dbmag[dgood]
dvmag= dvmag[dgood]
dvy  = dvy  [dgood]
dvx  = dvx  [dgood]
dvz  = dvz  [dgood]
dt   = dt   [dgood]
dden = dden [dgood]
dx   = dx   [dgood]
dy   = dy   [dgood]
dz   = dz   [dgood]

;Remove bad information WIND
wgood = where((wbx gt fillval) and (wby gt fillval) and (wbz gt fillval) and (wbmag gt fillval) and (wvmag gt fillval) $
    and (wvx gt fillval) and (wvy gt fillval) and (wvz gt fillval) and (wt gt fillval) and (wden gt fillval) $
    and (wx gt fillval) and (wy gt fillval) and (wz gt fillval))
wdoy = wdoy [wgood]
wbx  = wbx  [wgood]
wby  = wby  [wgood]
wbz  = wbz  [wgood]
wbmag= wbmag[wgood]
wvmag= wvmag[wgood]
wvy  = wvy  [wgood]
wvx  = wvx  [wgood]
wvz  = wvz  [wgood]
wt   = wt   [wgood]
wden = wden [wgood]
wx   = wx   [wgood]
wy   = wy   [wgood]
wz   = wz   [wgood]

;Remove bad information ACE
agood = where((abx gt fillval) and (aby gt fillval) and (abz gt fillval) and (abmag gt fillval) and (avmag gt fillval) $
    and (avx gt fillval) and (avy gt fillval) and (avz gt fillval) and (at gt fillval) and (aden gt fillval) $
    and (ax gt fillval) and (ay gt fillval) and (az gt fillval))
adoy = adoy [agood]
abx  = abx  [agood]
aby  = aby  [agood]
abz  = abz  [agood]
abmag= abmag[agood]
avmag= avmag[agood]
avy  = avy  [agood]
avx  = avx  [agood]
avz  = avz  [agood]
at   = at   [agood]
aden = aden [agood]
ax   = ax   [agood]
ay   = ay   [agood]
az   = az   [agood]

;Remove bad information DSCOVR
igood = where((ibx gt fillval) and (iby gt fillval) and (ibz gt fillval) and (ibmag gt fillval) and (ivmag gt fillval) $
    and (ivx gt fillval) and (ivy gt fillval) and (ivz gt fillval) and (it gt fillval) and (iden gt fillval) $
    and (ix gt fillval) and (iy gt fillval) and (iz gt fillval))
idoy = idoy [igood]
ibx  = ibx  [igood]
iby  = iby  [igood]
ibz  = ibz  [igood]
ibmag= ibmag[igood]
ivmag= ivmag[igood]
ivy  = ivy  [igood]
ivx  = ivx  [igood]
ivz  = ivz  [igood]
it   = it   [igood]
iden = iden [igood]
ix   = ix   [igood]
iy   = iy   [igood]
iz   = iz   [igood]




;interpolate wind onto dscovr
wbx   = interpol(temporary(wbx),wdoy,ddoy)
wby   = interpol(temporary(wby),wdoy,ddoy)
wbz   = interpol(temporary(wbz),wdoy,ddoy)
wvx   = interpol(temporary(wvx),wdoy,ddoy)
wvy   = interpol(temporary(wvy),wdoy,ddoy)
wvz   = interpol(temporary(wvz),wdoy,ddoy)
wx    = interpol(temporary(wx),wdoy,ddoy)
wy    = interpol(temporary(wy),wdoy,ddoy)
wz    = interpol(temporary(wz),wdoy,ddoy)
wbmag = interpol(temporary(wbmag),wdoy,ddoy)
wvmag = interpol(temporary(wvmag),wdoy,ddoy)
wt    = interpol(temporary(wt),wdoy,ddoy)
wden  = interpol(temporary(wden),wdoy,ddoy)

;interpolate wind onto dscovr
ibx   = interpol(temporary(ibx)  ,idoy,adoy)
iby   = interpol(temporary(iby)  ,idoy,adoy)
ibz   = interpol(temporary(ibz)  ,idoy,adoy)
ivx   = interpol(temporary(ivx)  ,idoy,adoy)
ivy   = interpol(temporary(ivy)  ,idoy,adoy)
ivz   = interpol(temporary(ivz)  ,idoy,adoy)
ix    = interpol(temporary(ix)   ,idoy,adoy)
iy    = interpol(temporary(iy)   ,idoy,adoy)
iz    = interpol(temporary(iz)   ,idoy,adoy)
ibmag = interpol(temporary(ibmag),idoy,adoy)
ivmag = interpol(temporary(ivmag),idoy,adoy)
it    = interpol(temporary(it)   ,idoy,adoy)
iden  = interpol(temporary(iden) ,idoy,adoy)


;Differences in parameters dscovr-wind
ddwbx   = dbx-wbx
ddwby   = dby-wby
ddwbz   = dbz-wbz
ddwvx   = dvx-wvx
ddwvy   = dvy-wvy
ddwvz   = dvz-wvz
ddwx    = dx-wx
ddwy    = dy-wy
ddwz    = dz-wz
ddwbmag = dbmag-wbmag
ddwvmag = dvmag-wvmag
ddwt    = dt-wt
ddwden  = dden-wden

;Differences in parameters ace-wind
daibx   = abx-ibx
daiby   = aby-iby
daibz   = abz-ibz
daivx   = avx-ivx
daivy   = avy-ivy
daivz   = avz-ivz
daix    = ax-ix
daiy    = ay-iy
daiz    = az-iz
daibmag = abmag-ibmag
daivmag = avmag-ivmag
dait    = at-it
daiden  = aden-iden


;standard deviation in dscovr-wind
sigddwbx   = stddev(ddwbx  )
sigddwby   = stddev(ddwby  )
sigddwbz   = stddev(ddwbz  )
sigddwvx   = stddev(ddwvx  )
sigddwvy   = stddev(ddwvy  )
sigddwvz   = stddev(ddwvz  )
sigddwx    = stddev(ddwx   )
sigddwy    = stddev(ddwy   )
sigddwz    = stddev(ddwz   )
sigddwbmag = stddev(ddwbmag)
sigddwvmag = stddev(ddwvmag)
sigddwt    = stddev(ddwt   )
sigddwden  = stddev(ddwden )

;standard deviation in ace-wind
sigdaibx   = stddev(daibx  ) 
sigdaiby   = stddev(daiby  ) 
sigdaibz   = stddev(daibz  ) 
sigdaivx   = stddev(daivx  ) 
sigdaivy   = stddev(daivy  ) 
sigdaivz   = stddev(daivz  ) 
sigdaix    = stddev(daix   ) 
sigdaiy    = stddev(daiy   ) 
sigdaiz    = stddev(daiz   ) 
sigdaibmag = stddev(daibmag) 
sigdaivmag = stddev(daivmag) 
sigdait    = stddev(dait   ) 
sigdaiden  = stddev(daiden ) 


;sigma
sigma = "162B ;" prevent byte from visually messing with color coding
plmns = "260B ;"
mu    = "154B ;"

sigma = "!4"+string(sigma)+"!X"
plmns = "!4"+string(plmns)+"!X"
mu    = "!4"+string(mu   )+"!X"

sigfmt = '("=",F5.1)'

print,mean(ddwvx)
;string for DSCOVR sigmas
strsigddwbx   =  mu+string([ mean(ddwbx  )],format=sigfmt)+", med"+string([ median(ddwbx  )],format=sigfmt)+", "+sigma+string([ sigddwbx  ],format=sigfmt)
strsigddwby   =  mu+string([ mean(ddwby  )],format=sigfmt)+", med"+string([ median(ddwby  )],format=sigfmt)+", "+sigma+string([ sigddwby  ],format=sigfmt)
strsigddwbz   =  mu+string([ mean(ddwbz  )],format=sigfmt)+", med"+string([ median(ddwbz  )],format=sigfmt)+", "+sigma+string([ sigddwbz  ],format=sigfmt)
strsigddwvx   =  mu+string([ mean(ddwvx  )],format=sigfmt)+", med"+string([ median(ddwvx  )],format=sigfmt)+", "+sigma+string([ sigddwvx  ],format=sigfmt)
strsigddwvy   =  mu+string([ mean(ddwvy  )],format=sigfmt)+", med"+string([ median(ddwvy  )],format=sigfmt)+", "+sigma+string([ sigddwvy  ],format=sigfmt)
strsigddwvz   =  mu+string([ mean(ddwvz  )],format=sigfmt)+", med"+string([ median(ddwvz  )],format=sigfmt)+", "+sigma+string([ sigddwvz  ],format=sigfmt)
strsigddwx    =  mu+string([ mean(ddwx   )],format=sigfmt)+", med"+string([ median(ddwx   )],format=sigfmt)+", "+sigma+string([ sigddwx   ],format=sigfmt)
strsigddwy    =  mu+string([ mean(ddwy   )],format=sigfmt)+", med"+string([ median(ddwy   )],format=sigfmt)+", "+sigma+string([ sigddwy   ],format=sigfmt)
strsigddwz    =  mu+string([ mean(ddwz   )],format=sigfmt)+", med"+string([ median(ddwz   )],format=sigfmt)+", "+sigma+string([ sigddwz   ],format=sigfmt)
strsigddwbmag =  mu+string([ mean(ddwbmag)],format=sigfmt)+", med"+string([ median(ddwbmag)],format=sigfmt)+", "+sigma+string([ sigddwbmag],format=sigfmt)
strsigddwvmag =  mu+string([ mean(ddwvmag)],format=sigfmt)+", med"+string([ median(ddwvmag)],format=sigfmt)+", "+sigma+string([ sigddwvmag],format=sigfmt)
strsigddwt    =  mu+string([ mean(ddwt   )],format=sigfmt)+", med"+string([ median(ddwt   )],format=sigfmt)+", "+sigma+string([ sigddwt   ],format=sigfmt)
strsigddwden  =  mu+string([ mean(ddwden )],format=sigfmt)+", med"+string([ median(ddwden )],format=sigfmt)+", "+sigma+string([ sigddwden ],format=sigfmt)
                                                                                                        
;string for ACE  sigmas                                       
strsigdaibx   =  mu+string([ mean(daibx  )],format=sigfmt)+", med"+string([ median(daibx  )],format=sigfmt)+", "+sigma+string([ sigdaibx  ],format=sigfmt)
strsigdaiby   =  mu+string([ mean(daiby  )],format=sigfmt)+", med"+string([ median(daiby  )],format=sigfmt)+", "+sigma+string([ sigdaiby  ],format=sigfmt)
strsigdaibz   =  mu+string([ mean(daibz  )],format=sigfmt)+", med"+string([ median(daibz  )],format=sigfmt)+", "+sigma+string([ sigdaibz  ],format=sigfmt)
strsigdaivx   =  mu+string([ mean(daivx  )],format=sigfmt)+", med"+string([ median(daivx  )],format=sigfmt)+", "+sigma+string([ sigdaivx  ],format=sigfmt)
strsigdaivy   =  mu+string([ mean(daivy  )],format=sigfmt)+", med"+string([ median(daivy  )],format=sigfmt)+", "+sigma+string([ sigdaivy  ],format=sigfmt)
strsigdaivz   =  mu+string([ mean(daivz  )],format=sigfmt)+", med"+string([ median(daivz  )],format=sigfmt)+", "+sigma+string([ sigdaivz  ],format=sigfmt)
strsigdaix    =  mu+string([ mean(daix   )],format=sigfmt)+", med"+string([ median(daix   )],format=sigfmt)+", "+sigma+string([ sigdaix   ],format=sigfmt)
strsigdaiy    =  mu+string([ mean(daiy   )],format=sigfmt)+", med"+string([ median(daiy   )],format=sigfmt)+", "+sigma+string([ sigdaiy   ],format=sigfmt)
strsigdaiz    =  mu+string([ mean(daiz   )],format=sigfmt)+", med"+string([ median(daiz   )],format=sigfmt)+", "+sigma+string([ sigdaiz   ],format=sigfmt)
strsigdaibmag =  mu+string([ mean(daibmag)],format=sigfmt)+", med"+string([ median(daibmag)],format=sigfmt)+", "+sigma+string([ sigdaibmag],format=sigfmt)
strsigdaivmag =  mu+string([ mean(daivmag)],format=sigfmt)+", med"+string([ median(daivmag)],format=sigfmt)+", "+sigma+string([ sigdaivmag],format=sigfmt)
strsigdait    =  mu+string([ mean(dait   )],format=sigfmt)+", med"+string([ median(dait   )],format=sigfmt)+", "+sigma+string([ sigdait   ],format=sigfmt)
strsigdaiden  =  mu+string([ mean(daiden )],format=sigfmt)+", med"+string([ median(daiden )],format=sigfmt)+", "+sigma+string([ sigdaiden ],format=sigfmt)







;number of samples for normalization
samples = float(n_elements(wdoy))
;ace samples
aamples = float(n_elements(adoy))


;Set up plots
device,decomposed=0
device,retain=2
loadct,12
DLM_LOAD,'PNG'



dlet = "104B;" prevent byte from visually messing with color coding 
delta = '!4'+string(dlet)+'!X'
!p.thick=5
!x.thick=4
!y.thick=4
!x.charsize=2.2
!y.charsize=2.2
!p.charthick=2.2
!p.color=0
!p.background=255
!y.style = 1
!y.minor = 2
!y.ticklen = .04
!x.ticklen = .04

;Set up histograms
vxsize = 10
vysize = 10
vzsize = 10
desize = 0.5
thsize = 5
vxrange = [-100,100]
vyrange = [0,50]

;locations for xyouts
xxylocd = -95 
yxylocd = 45
yxyloca = 42


;Measured Velocity components dscovr
vxdwhist = histogram(ddwvx,binsize=vxsize,location=vxdwbins,min=-100,max=100)
format_hist,vxdwhist,vxdwbins,vxsize
vydwhist = histogram(ddwvy,binsize=vysize,location=vydwbins,min=-100,max=100)
format_hist,vydwhist,vydwbins,vysize
vzdwhist = histogram(ddwvz,binsize=vzsize,location=vzdwbins,min=-100,max=100)
format_hist,vzdwhist,vzdwbins,vzsize

;Measured Velocity components  ace
vxaihist = histogram(daivx,binsize=vxsize,location=vxaibins,min=-100,max=100)
format_hist,vxaihist,vxaibins,vxsize
vyaihist = histogram(daivy,binsize=vysize,location=vyaibins,min=-100,max=100)
format_hist,vyaihist,vyaibins,vysize
vzaihist = histogram(daivz,binsize=vzsize,location=vzaibins,min=-100,max=100)
format_hist,vzaihist,vzaibins,vzsize

;plot measured velocity components
plot,vxdwbins,vxdwhist,psym=10, $
    xtitle=delta+'Vx (X-WIND) [km/s]',ytitle='Occurence [%]',position=[.08,.6,.32,.95], $
    xrange=vxrange,yrange=vyrange
    oplot,vxaibins,vxaihist,psym=10,color=200
    xyouts,xxylocd,yxylocd,strsigddwvx,charsize=1.5,charthick=2.5
    xyouts,xxylocd,yxyloca,strsigdaivx,charsize=1.5,charthick=2.5,color=200
    xyouts,55,yxylocd,'(DSCOVR)',charsize=1.5,charthick=2.5
    xyouts,55,yxyloca,'(ACE)',charsize=1.5,charthick=2.5,color=200

plot,vydwbins,vydwhist,psym=10,/NOERASE, $
    xtitle=delta+'Vy (X-WIND) [km/s]',ytitle='Occurence [%]',position=[.41,.6,.65,.95], $
    xrange=vxrange,yrange=vyrange
    oplot,vyaibins,100.*vyaihist/float(total(vyaihist)),psym=10,color=200
    xyouts,xxylocd,yxylocd,strsigddwvy,charsize=1.5,charthick=2.5
    xyouts,xxylocd,yxyloca,strsigdaivy,charsize=1.5,charthick=2.5,color=200

plot,vzdwbins,vzdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Vz (X-WIND) [km/s]',ytitle='Occurence [%]',position=[.74,.6,.98,.95], $
    xrange=vxrange,yrange=vyrange
    oplot,vzaibins,vzaihist,psym=10,color=200
    xyouts,xxylocd,yxylocd,strsigddwvz,charsize=1.5,charthick=2.5
    xyouts,xxylocd,yxyloca,strsigdaivz,charsize=1.5,charthick=2.5,color=200


;Measured Velocity components 
vmagdwhist = histogram(ddwvmag,binsize=vxsize,location=vmagdwbins,min=-100,max=100)
format_hist,vmagdwhist,vmagdwbins,vxsize
tdwhist = histogram(ddwt,binsize=thsize,location=tdwbins,min=-100,max=100)
format_hist,tdwhist,tdwbins,thsize
dendwhist = histogram(ddwden,binsize=desize,location=dendwbins,min=-10,max=10)
format_hist,dendwhist,dendwbins,desize

;Measured Velocity components 
vmagaihist = histogram(daivmag,binsize=vxsize,location=vmagaibins,min=-100,max=100)
format_hist,vmagaihist,vmagaibins,vxsize
taihist = histogram(dait,binsize=thsize,location=taibins,min=-100,max=100)
format_hist,taihist,taibins,thsize
denaihist = histogram(daiden,binsize=desize,location=denaibins,min=-10,max=10)
format_hist,denaihist,denaibins,desize

;plot measured velocity components
plot,vmagdwbins,vmagdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Speed (X-WIND) [km/s]',ytitle='Occurence [%]',position=[.08,.10,.32,.45], $
    xrange=vxrange,yrange=vyrange
    oplot,vmagaibins,vmagaihist,psym=10,color=200
    xyouts,xxylocd,yxylocd,strsigddwvmag,charsize=1.5,charthick=2.5
    xyouts,xxylocd,yxyloca,strsigdaivmag,charsize=1.5,charthick=2.5,color=200

plot,tdwbins,tdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Th. Speed (X-WIND) [km/s]',ytitle='Occurence [%]',position=[.41,.10,.65,.45], $
    xrange=[-50,50],yrange=vyrange
    oplot,taibins,taihist,psym=10,color=200
    xyouts,-55,yxylocd,strsigddwt,charsize=1.5,charthick=2.5
    xyouts,-55,yxyloca,strsigdait,charsize=1.5,charthick=2.5,color=200

plot,dendwbins,dendwhist,psym=10,/NOERASE, $
    xtitle=delta+'Den. (X-WIND) [cm^-3]',ytitle='Occurence [%]',position=[.74,.10,.98,.45], $
    xrange=[-5,5],yrange=vyrange
    oplot,denaibins,denaihist,psym=10,color=200
    xyouts,-5.2,yxylocd,strsigddwden,charsize=1.5,charthick=2.5
    xyouts,-5.2,yxyloca,strsigdaiden,charsize=1.5,charthick=2.5,color=200



write_png,'../compare_plots/compare_grid.png',tvrd(/true)


end
