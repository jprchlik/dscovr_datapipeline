;----------------------------------------------------------------------------------------
;
;Uses
;
;USAGE 
;compare_wind_dscovr_ace,syear,sdoy,edoy=edoy,aceoff=aceoff,eyear=eyear
;
;
;----------------------------------------------------------------------------------------
FUNCTION Exponent, axis, index, number

  ; A special case.
  IF number EQ 0 THEN RETURN, '0' 

  ; Assuming multiples of 10 with format.
  ex = String(number, Format='(e8.0)') 
  pt = StrPos(ex, '.')

  first = StrMid(ex, 0, pt)
  sign = StrMid(ex, pt+2, 1)
  thisExponent = StrMid(ex, pt+3)

  ; Shave off leading zero in exponent
  WHILE StrMid(thisExponent, 0, 1) EQ '0' DO thisExponent = StrMid(thisExponent, 1)

  ; Fix for sign and missing zero problem.
  IF (Long(thisExponent) EQ 0) THEN BEGIN
     sign = ''
     thisExponent = '0'
  ENDIF

  ; Make the exponent a superscript.
  IF sign EQ '-' THEN BEGIN
     RETURN, '10!U' + sign + thisExponent + '!N' 
  ENDIF ELSE BEGIN             
     RETURN, '10!U' + thisExponent + '!N'
  ENDELSE
  
END



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
if keyword_set(span) then span = span else span = 25 ;span to loop +/- in minutes
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
    offset = ddoy+tvals[i]
    sdx = interpol(wx,wdoy,offset)
    sdx = temporary(sdx)-dx
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
chim = where(fcvals eq min(fcvals,/nan))
if n_elements(size(chim)) gt 3 then begin ;move the time as little as possible
    dtime = ftvals[chim]
    minmo = where(abs(dtime) eq min(abs(dtime),/nan)) ; move the minimum amount
    dtime = dtime[minmo]
endif else begin 
    dtime = 0. ; return 0 if no min??
endelse
;print,'Chi^2 min value',strcompress(min(fcvals),/remove_all)
;print,'Chi^2 time offset',strcompress(dtime*24.*60.),'min'


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
    rep = where(iday eq days[i])
    deldoy = chi_min_time(wdoy[wcut],wvx[wcut],doy[cut],vx[cut],ddx=dvx[cut])
    new_doy[rep] = deldoy;add time offset to new array
endfor

doy = doy+new_doy 

return,doy
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



;----------------------------------------------------------------------------------------
;USAGE
;running_med,x,y,medx,medy,bins=bins,sig=sig
;----------------------------------------------------------------------------------------
pro running_med,x,y,medx,medy,bins=bins,sig=sig
if keyword_set(bins) then bins=bins else bins = 10

;Find the number of bins to create the running median
maxx = max(x)
minx = 0.
ranx = maxx-minx
binx = ranx/bins
numx = fix(binx)+1

;Create arrays which will contain the median values
medy = fltarr(numx)
medx = fltarr(numx)
sig  = fltarr(numx)

;loop over all bin values
for i=0,numx-1 do begin
    medx[i] = float(i)*float(bins)+float(bins)/2.
    good_val = where((x ge float(i)*float(bins)) and (x lt float(i+1.)*float(bins)),count)
    if count gt 0 then begin
        medy[i] = median(y[good_val])
        sig [i] = stddev(y[good_val])/float(n_elements(good_val))
    endif else begin
        medy[i] = -9999.0
        sig [i] = -9999.0
    endelse
endfor
end

;--------------------------------------------------
;USAGE
;    error_bars,xval,yval,err,p_err_x,p_err_y
;
;COMMENTS
;    creates error bars for plots 
;    Allow the syntax claims only Y errors you can flop x and y to plot x errors
;--------------------------------------------------
pro error_bars,xvals,yvals,svals,ecol=ecol
if keyword_set(ecol) then ecol=ecol else ecol=200
for i=0,n_elements(xvals)-1 do begin
    xval = xvals[i]
    yval = yvals[i]
    sval = svals[i]
    ex = [xval,xval]
    ey = [abs(sval),-abs(sval)]+yval
    if yval gt -9998. then oplot,ex,ey,color=ecol,linestyle=0,thick=.2
endfor
end


;----------------------------------------------------------------------------------------
;MAIN LOOP
;----------------------------------------------------------------------------------------
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
        ddoy = time_offset_cor(wdoy,wvx,temporary(ddoy),dvx)

;do the same for ace
        load_plasma,ayear[i],tsdoy,tedoy,/ACE,doy=adoy,bx=abx,by=aby,bz=abz,bmag=abmag,vmag=avmag,vy=avy,vx=avx,vz=avz,t=at,den=aden,x=ax,y=ay,z=az
        load_plasma,ayear[i],tsdoy,tedoy,/WIND,doy=idoy,bx=ibx,by=iby,bz=ibz,bmag=ibmag,vmag=ivmag,vy=ivy,vx=ivx,vz=ivz,t=it,den=iden,x=ix,y=iy,z=iz
;correct for orbital time offsets ACE/WIND
        adoy = time_offset_cor(idoy,ivx,temporary(adoy),avx)

    endif else begin
         tsdoy = 1
         load_plasma,ryear[i],tsdoy,tedoy,/DSC,doy=tddoy,bx=tdbx,by=tdby,bz=tdbz,bmag=tdbmag,vmag=tdvmag,vy=tdvy,vx=tdvx,vz=tdvz,t=tdt,den=tdden,x=tdx,y=tdy,z=tdz
         load_plasma,ryear[i],tsdoy,tedoy,/WIND,doy=twdoy,bx=twbx,by=twby,bz=twbz,bmag=twbmag,vmag=twvmag,vy=twvy,vx=twvx,vz=twvz,t=twt,den=twden,x=twx,y=twy,z=twz
;correct for orbital time offsets DSCOVR/WIND
         print,'DSCOVR 2'
         tddoy = time_offset_cor(twdoy,twvx,temporary(tddoy),tdvx)

;do the same for ace
         load_plasma,ayear[i],tsdoy,tedoy,/ACE,doy=tadoy,bx=tabx,by=taby,bz=tabz,bmag=tabmag,vmag=tavmag,vy=tavy,vx=tavx,vz=tavz,t=tat,den=taden,x=tax,y=tay,z=taz
         load_plasma,ayear[i],tsdoy,tedoy,/WIND,doy=tidoy,bx=tibx,by=tiby,bz=tibz,bmag=tibmag,vmag=tivmag,vy=tivy,vx=tivx,vz=tivz,t=tit,den=tiden,x=tix,y=tiy,z=tiz
;correct for orbital time offsets ACE/WIND
         print,'ACE 2'
         tadoy = time_offset_cor(tidoy,tivx,temporary(tadoy),tavx)



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
dgood = where((dvmag gt fillval) $
    and (dvx gt fillval) and (dvy gt fillval) and (dvz gt fillval) and (dt gt fillval) and (dden gt fillval) and (dden le 1000) ) ;max to kill density value of 1.9E10
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
wgood = where((wvmag gt fillval) $
    and (wvx gt fillval) and (wvy gt fillval) and (wvz gt fillval) and (wt gt fillval) and (wden gt fillval)) 
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
agood = where((avmag gt fillval) $
    and (avx gt fillval) and (avy gt fillval) and (avz gt fillval) and (at gt fillval) and (aden gt fillval) )
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

;Remove bad information WIND
igood = where( (ivmag gt fillval) $
    and (ivx gt fillval) and (ivy gt fillval) and (ivz gt fillval) and (it gt fillval) and (iden gt fillval)) 
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
ddwt    = 100.*(dt-wt)/wt
ddwden  = 100.*(dden-wden)/wden

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
dait    = 100.*(at-it)/it
daiden  = 100.*(aden-iden)/iden


;standard deviation in dscovr-wind
;Cut a 99.7 percentile (3sigma)
perarray = [0.003,0.997]

;Find percentile values
perddwbx   = cgpercentiles(ddwbx  ,percentiles=perarray)
perddwby   = cgpercentiles(ddwby  ,percentiles=perarray)
perddwbz   = cgpercentiles(ddwbz  ,percentiles=perarray)
perddwvx   = cgpercentiles(ddwvx  ,percentiles=perarray)
perddwvy   = cgpercentiles(ddwvy  ,percentiles=perarray)
perddwvz   = cgpercentiles(ddwvz  ,percentiles=perarray)
perddwx    = cgpercentiles(ddwx   ,percentiles=perarray)
perddwy    = cgpercentiles(ddwy   ,percentiles=perarray)
perddwz    = cgpercentiles(ddwz   ,percentiles=perarray)
perddwbmag = cgpercentiles(ddwbmag,percentiles=perarray)
perddwvmag = cgpercentiles(ddwvmag,percentiles=perarray)
perddwt    = cgpercentiles(ddwt   ,percentiles=perarray)
perddwden  = cgpercentiles(ddwden ,percentiles=perarray)

;create array slices
cutddwbx   = where((ddwbx   gt perddwbx   [0]) and (ddwbx   lt perddwbx   [1]))
cutddwby   = where((ddwby   gt perddwby   [0]) and (ddwby   lt perddwby   [1]))
cutddwbz   = where((ddwbz   gt perddwbz   [0]) and (ddwbz   lt perddwbz   [1]))
cutddwvx   = where((ddwvx   gt perddwvx   [0]) and (ddwvx   lt perddwvx   [1]))
cutddwvy   = where((ddwvy   gt perddwvy   [0]) and (ddwvy   lt perddwvy   [1]))
cutddwvz   = where((ddwvz   gt perddwvz   [0]) and (ddwvz   lt perddwvz   [1]))
cutddwx    = where((ddwx    gt perddwx    [0]) and (ddwx    lt perddwx    [1]))
cutddwy    = where((ddwy    gt perddwy    [0]) and (ddwy    lt perddwy    [1]))
cutddwz    = where((ddwz    gt perddwz    [0]) and (ddwz    lt perddwz    [1]))
cutddwbmag = where((ddwbmag gt perddwbmag [0]) and (ddwbmag lt perddwbmag [1]))
cutddwvmag = where((ddwvmag gt perddwvmag [0]) and (ddwvmag lt perddwvmag [1]))
cutddwt    = where((ddwt    gt perddwt    [0]) and (ddwt    lt perddwt    [1]))
cutddwden  = where((ddwden  gt perddwden  [0]) and (ddwden  lt perddwden  [1]))


;calculate cut standard deviation
sigddwbx   = stddev(ddwbx  [cutddwbx   ])
sigddwby   = stddev(ddwby  [cutddwby   ])
sigddwbz   = stddev(ddwbz  [cutddwbz   ])
sigddwvx   = stddev(ddwvx  [cutddwvx   ])
sigddwvy   = stddev(ddwvy  [cutddwvy   ])
sigddwvz   = stddev(ddwvz  [cutddwvz   ])
sigddwx    = stddev(ddwx   [cutddwx    ])
sigddwy    = stddev(ddwy   [cutddwy    ])
sigddwz    = stddev(ddwz   [cutddwz    ])
sigddwbmag = stddev(ddwbmag[cutddwbmag ])
sigddwvmag = stddev(ddwvmag[cutddwvmag ])
sigddwt    = stddev(ddwt   [cutddwt    ])
sigddwden  = stddev(ddwden [cutddwden  ])

;standard deviation in ace-wind
;Find percentile values
perdaibx   = cgpercentiles(daibx  ,percentiles=perarray)
perdaiby   = cgpercentiles(daiby  ,percentiles=perarray)
perdaibz   = cgpercentiles(daibz  ,percentiles=perarray)
perdaivx   = cgpercentiles(daivx  ,percentiles=perarray)
perdaivy   = cgpercentiles(daivy  ,percentiles=perarray)
perdaivz   = cgpercentiles(daivz  ,percentiles=perarray)
perdaix    = cgpercentiles(daix   ,percentiles=perarray)
perdaiy    = cgpercentiles(daiy   ,percentiles=perarray)
perdaiz    = cgpercentiles(daiz   ,percentiles=perarray)
perdaibmag = cgpercentiles(daibmag,percentiles=perarray)
perdaivmag = cgpercentiles(daivmag,percentiles=perarray)
perdait    = cgpercentiles(dait   ,percentiles=perarray)
perdaiden  = cgpercentiles(daiden ,percentiles=perarray)

;create array slices
cutdaibx   = where((daibx   gt perdaibx   [0]) and (daibx   lt perdaibx   [1]))
cutdaiby   = where((daiby   gt perdaiby   [0]) and (daiby   lt perdaiby   [1]))
cutdaibz   = where((daibz   gt perdaibz   [0]) and (daibz   lt perdaibz   [1]))
cutdaivx   = where((daivx   gt perdaivx   [0]) and (daivx   lt perdaivx   [1]))
cutdaivy   = where((daivy   gt perdaivy   [0]) and (daivy   lt perdaivy   [1]))
cutdaivz   = where((daivz   gt perdaivz   [0]) and (daivz   lt perdaivz   [1]))
cutdaix    = where((daix    gt perdaix    [0]) and (daix    lt perdaix    [1]))
cutdaiy    = where((daiy    gt perdaiy    [0]) and (daiy    lt perdaiy    [1]))
cutdaiz    = where((daiz    gt perdaiz    [0]) and (daiz    lt perdaiz    [1]))
cutdaibmag = where((daibmag gt perdaibmag [0]) and (daibmag lt perdaibmag [1]))
cutdaivmag = where((daivmag gt perdaivmag [0]) and (daivmag lt perdaivmag [1]))
cutdait    = where((dait    gt perdait    [0]) and (dait    lt perdait    [1]))
cutdaiden  = where((daiden  gt perdaiden  [0]) and (daiden  lt perdaiden  [1]))


;calculate cut standard deviation
sigdaibx   = stddev(daibx  [cutdaibx   ])
sigdaiby   = stddev(daiby  [cutdaiby   ])
sigdaibz   = stddev(daibz  [cutdaibz   ])
sigdaivx   = stddev(daivx  [cutdaivx   ])
sigdaivy   = stddev(daivy  [cutdaivy   ])
sigdaivz   = stddev(daivz  [cutdaivz   ])
sigdaix    = stddev(daix   [cutdaix    ])
sigdaiy    = stddev(daiy   [cutdaiy    ])
sigdaiz    = stddev(daiz   [cutdaiz    ])
sigdaibmag = stddev(daibmag[cutdaibmag ])
sigdaivmag = stddev(daivmag[cutdaivmag ])
sigdait    = stddev(dait   [cutdait    ])
sigdaiden  = stddev(daiden [cutdaiden  ])

;standard deviation in ace-wind
;sigdaibx   = stddev(daibx  ) 
;sigdaiby   = stddev(daiby  ) 
;sigdaibz   = stddev(daibz  ) 
;sigdaibx   = stddev(daibx  ) 
;sigdaiby   = stddev(daiby  ) 
;sigdaibz   = stddev(daibz  ) 
;sigdaivx   = stddev(daivx  ) 
;sigdaivy   = stddev(daivy  ) 
;sigdaivz   = stddev(daivz  ) 
;sigdaix    = stddev(daix   ) 
;sigdaiy    = stddev(daiy   ) 
;sigdaiz    = stddev(daiz   ) 
;sigdaibmag = stddev(daibmag) 
;sigdaivmag = stddev(daivmag) 
;sigdait    = stddev(dait   ) 
;sigdaiden  = stddev(daiden ) 


;sigma
;sigma = "162B ;" prevent byte from visually messing with color coding
;plmns = "260B ;"
;mu    = "154B ;"
;
;sigma = "!4"+string(sigma)+"!X"
;plmns = "!4"+string(plmns)+"!X"
;mu    = "!4"+string(mu   )+"!X"

sigma = cgsymbol('sigma')
plmns = cgsymbol('+-')
mu    = cgsymbol('mu')

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
print,mean(ddwden),sigddwden
                                                                                                        
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
set_plot,'ps'
;device,decomposed=0
;device,retain=2
loadct,12
;DLM_LOAD,'PNG'



;dlet = "104B;" prevent byte from visually messing with color coding 
delta = cgsymbol('delta');'!4'+string(dlet)+'!X'
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
!P.FONT = 0

;Set up histograms
vxsize = 10
vysize = 10
vzsize = 10
desize = 5
thsize = 5
vxrange = [-100,100]
vyrange = [0,55]

;locations for xyouts
xxylocd = -95 
yxylocd = 47
yxyloca = 43

;plot locations
plot1 = [.08,.6,.32,.95]
plot2 = [.41,.6,.65,.95]
plot3 = [.74,.6,.98,.95]
plot4 = [.08,.10,.32,.45]
plot5 = [.41,.10,.65,.45] 
plot6 = [.74,.10,.98,.45]

;-----------------------------------------------------------------------------
;Compare distributions for ACE and DSCOVR
;-----------------------------------------------------------------------------
device,filename='../compare_plots/compare_grid.eps',encapsulated=1,/helvetica,xsize=8.0,ysize=5.3,/inch

;Measured Velocity components dscovr
vxdwhist = histogram(ddwvx,binsize=vxsize,location=vxdwbins,min=-1000,max=1000)
format_hist,vxdwhist,vxdwbins,vxsize
vydwhist = histogram(ddwvy,binsize=vysize,location=vydwbins,min=-1000,max=1000)
format_hist,vydwhist,vydwbins,vysize
vzdwhist = histogram(ddwvz,binsize=vzsize,location=vzdwbins,min=-1000,max=1000)
format_hist,vzdwhist,vzdwbins,vzsize

;Measured Velocity components  ace
vxaihist = histogram(daivx,binsize=vxsize,location=vxaibins,min=-1000,max=1000)
format_hist,vxaihist,vxaibins,vxsize
vyaihist = histogram(daivy,binsize=vysize,location=vyaibins,min=-1000,max=1000)
format_hist,vyaihist,vyaibins,vysize
vzaihist = histogram(daivz,binsize=vzsize,location=vzaibins,min=-1000,max=1000)
format_hist,vzaihist,vzaibins,vzsize

;plot measured velocity components
plot,vxdwbins,vxdwhist,psym=10, $
    xtitle=delta+'Vx [km/s]',ytitle='Occurence [%]',position=plot1, $
    xrange=vxrange,yrange=vyrange
    oplot,vxaibins,vxaihist,psym=10,color=200
    xyouts,xxylocd,yxylocd,strsigddwvx,charsize=0.7,charthick=2.5
    xyouts,xxylocd,yxyloca,strsigdaivx,charsize=0.7,charthick=2.5,color=200

plot,vydwbins,vydwhist,psym=10,/NOERASE, $
    xtitle=delta+'Vy [km/s]',ytitle='Occurence [%]',position=plot2, $
    xrange=vxrange,yrange=vyrange
    oplot,vyaibins,vyaihist,psym=10,color=200
    xyouts,xxylocd,yxylocd,strsigddwvy,charsize=0.7,charthick=2.5
    xyouts,xxylocd,yxyloca,strsigdaivy,charsize=0.7,charthick=2.5,color=200

plot,vzdwbins,vzdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Vz [km/s]',ytitle='Occurence [%]',position=plot3, $
    xrange=vxrange,yrange=vyrange
    oplot,vzaibins,vzaihist,psym=10,color=200
    xyouts,xxylocd,yxylocd,strsigddwvz,charsize=0.7,charthick=2.5
    xyouts,xxylocd,yxyloca,strsigdaivz,charsize=0.7,charthick=2.5,color=200


;Measured Velocity components 
vmagdwhist = histogram(ddwvmag,binsize=vxsize,location=vmagdwbins,min=-1000,max=1000)
format_hist,vmagdwhist,vmagdwbins,vxsize
tdwhist = histogram(ddwt,binsize=thsize,location=tdwbins,min=-1000,max=1000)
format_hist,tdwhist,tdwbins,thsize
dendwhist = histogram(ddwden,binsize=desize,location=dendwbins,min=-1000,max=1000)
format_hist,dendwhist,dendwbins,desize

;Measured Velocity components 
vmagaihist = histogram(daivmag,binsize=vxsize,location=vmagaibins,min=-1000,max=1000)
format_hist,vmagaihist,vmagaibins,vxsize
taihist = histogram(dait,binsize=thsize,location=taibins,min=-1000,max=1000)
format_hist,taihist,taibins,thsize
denaihist = histogram(daiden,binsize=desize,location=denaibins,min=-1000,max=1000)
format_hist,denaihist,denaibins,desize

;plot measured velocity components
plot,vmagdwbins,vmagdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Speed [km/s]',ytitle='Occurence [%]',position=plot4, $
    xrange=vxrange,yrange=vyrange
    oplot,vmagaibins,vmagaihist,psym=10,color=200
    xyouts,xxylocd,yxylocd,strsigddwvmag,charsize=0.7,charthick=2.5
    xyouts,xxylocd,yxyloca,strsigdaivmag,charsize=0.7,charthick=2.5,color=200

plot,tdwbins,tdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Th. Speed [%]',ytitle='Occurence [%]',position=plot5, $
    xrange=[-50,50],yrange=vyrange
    oplot,taibins,taihist,psym=10,color=200
    xyouts,-55,yxylocd,strsigddwt,charsize=0.7,charthick=2.5
    xyouts,-55,yxyloca,strsigdait,charsize=0.7,charthick=2.5,color=200

plot,dendwbins,dendwhist,psym=10,/NOERASE, $
    xtitle=delta+'Den. [%]',ytitle='Occurence [%]',position=plot6, $
    xrange=[-50,50],yrange=vyrange
    oplot,denaibins,denaihist,psym=10,color=200
    xyouts,-55.,yxylocd,strsigddwden,charsize=0.7,charthick=2.5
    xyouts,-55.,yxyloca,strsigdaiden,charsize=0.7,charthick=2.5,color=200

device,/close
;write_png,'../compare_plots/compare_grid.png',tvrd(/true)

;-----------------------------------------------------------------------------
;Compare log distributions for ACE and DSCOVR
;-----------------------------------------------------------------------------

;Set up eps plot
device,filename='../compare_plots/compare_grid_log.eps',encapsulated=1,/helvetica,xsize=8.5,ysize=5.0,/inch
!P.Charsize = .5
!P.charthick = 2.0
lvyrange = [1.e-4,100]

;yrange for 2 simga overplotting
ysiglim = [-10000,10000]
;plot measured velocity components
plot,vxdwbins,vxdwhist,psym=10, $
    xtitle=delta+'Vx [km/s]',ytitle='Occurence [%]',position=plot1, $
    xrange=[-600,600],yrange=lvyrange,xticks=4,xstyle=1,xminor=3, $
    ytickformat='exponent',ylog=1,ystyle=1,yminor=5

    oplot,vxaibins,vxaihist,psym=10,color=200
    oplot, sigddwvx*[2.,2.],ysiglim,psym=10,linestyle=2
    oplot,-sigddwvx*[2.,2.],ysiglim,psym=10,linestyle=2

plot,vydwbins,vydwhist,psym=10,/NOERASE, $
    xtitle=delta+'Vy [km/s]',ytitle='Occurence [%]',position=plot2, $
    xrange=vxrange*3,yrange=lvyrange, $
    ytickformat='exponent',ylog=1,ystyle=1,yminor=5
    oplot,vyaibins,vyaihist,psym=10,color=200
    oplot, sigddwvy*[2.,2.],ysiglim,psym=10,linestyle=2
    oplot,-sigddwvy*[2.,2.],ysiglim,psym=10,linestyle=2

plot,vzdwbins,vzdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Vz [km/s]',ytitle='Occurence [%]',position=plot3, $
    xrange=vxrange*3,yrange=lvyrange, $
    ytickformat='exponent',ylog=1,ystyle=1,yminor=5
    oplot,vzaibins,vzaihist,psym=10,color=200
    oplot, sigddwvz*[2.,2.],ysiglim,psym=10,linestyle=2
    oplot,-sigddwvz*[2.,2.],ysiglim,psym=10,linestyle=2


;plot measured velocity components
plot,vmagdwbins,vmagdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Speed [km/s]',ytitle='Occurence [%]',position=plot4, $
    xrange=[-600,600],yrange=lvyrange,xticks=4,xstyle=1,xminor=3, $
    ytickformat='exponent',ylog=1,ystyle=1,yminor=5
    oplot,vmagaibins,vmagaihist,psym=10,color=200
    oplot, sigddwvmag*[2.,2.],ysiglim,psym=10,linestyle=2
    oplot,-sigddwvmag*[2.,2.],ysiglim,psym=10,linestyle=2

plot,tdwbins,tdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Th. Speed [%]',ytitle='Occurence [%]',position=plot5, $
    xrange=[-200,600],yrange=lvyrange, $
    ytickformat='exponent',ylog=1,ystyle=1,yminor=5
    oplot,taibins,taihist,psym=10,color=200
    oplot, sigddwt*[2.,2.],ysiglim,psym=10,linestyle=2
    oplot,-sigddwt*[2.,2.],ysiglim,psym=10,linestyle=2

plot,dendwbins,dendwhist,psym=10,/NOERASE, $
    xtitle=delta+'Den. [%]',ytitle='Occurence [%]',position=plot6, $
    xrange=[-2,6]*100,yrange=lvyrange,xstyle=1, $
    ytickformat='exponent',ylog=1,ystyle=1,yminor=5
    oplot,denaibins,denaihist,psym=10,color=200
    oplot, sigddwden*[2.,2.],ysiglim,psym=10,linestyle=2
    oplot,-sigddwden*[2.,2.],ysiglim,psym=10,linestyle=2

legend,['DSC-WIND','ACE-WIND'],colors=[0,200],linestyle=[0,0],/right,box=0,charsize=0.65,charthick=2.0

device,/close
;write_png,'../compare_plots/compare_grid_log.png',tvrd(/true)

;-----------------------------------------------------------------------------
;Check that plot ranges are similar for ACE and DSCOVR
;-----------------------------------------------------------------------------
device,filename='../compare_plots/compare_phys_grid_log.eps',encapsulated=1,/helvetica
;write_png,'../compare_plots/compare_phys_grid.png',tvrd(/true)

;Measured Velocity components dscovr
vxdhist = histogram(dvx,binsize=20,location=vxdbins,min=-1000,max=1000)
format_hist,vxdhist,vxdbins,vxsize
vydhist = histogram(dvy,binsize=vysize,location=vydbins,min=-1000,max=1000)
format_hist,vydhist,vydbins,vysize
vzdhist = histogram(dvz,binsize=vzsize,location=vzdbins,min=-1000,max=1000)
format_hist,vzdhist,vzdbins,vzsize

;Measured Velocity components  ace
vxahist = histogram(avx,binsize=20,location=vxabins,min=-1000,max=1000)
format_hist,vxahist,vxabins,vxsize
vyahist = histogram(avy,binsize=vysize,location=vyabins,min=-1000,max=1000)
format_hist,vyahist,vyabins,vysize
vzahist = histogram(avz,binsize=vzsize,location=vzabins,min=-1000,max=1000)
format_hist,vzahist,vzabins,vzsize

;plot measured velocity components
yrange2 = [0.,25.]
plot,vxdbins,vxdhist,psym=10, $
    xtitle='Vx [km/s]',ytitle='Occurence [%]',position=plot1,$
    xrange=[-800.,-200.],yrange=yrange2,ystyle=1
    oplot,vxabins,vxahist,psym=10,color=200

plot,vydbins,vydhist,psym=10, /NOERASE, $
    xtitle='Vy [km/s]',ytitle='Occurence [%]',position=plot2,$
    xrange=[-200.,200.],yrange=yrange2,ystyle=1
    oplot,vyabins,vyahist,psym=10,color=200

plot,vzdbins,vzdhist,psym=10, /NOERASE, $
    xtitle='Vz [km/s]',ytitle='Occurence [%]',position=plot3,$
    xrange=[-200,200.],yrange=yrange2,ystyle=1
    oplot,vzabins,vzahist,psym=10,color=200

;Measured Physical components DSCOVR 
vmagdhist = histogram(dvmag,binsize=20,location=vmagdbins,min=-1000,max=1000)
format_hist,vmagdhist,vmagdbins,vxsize
tdhist = histogram(dt,binsize=thsize,location=tdbins,min=-1000,max=1000)
format_hist,tdhist,tdbins,thsize
dendhist = histogram(dden,binsize=desize,location=dendbins,min=-1000,max=1000)
format_hist,dendhist,dendbins,desize

;Measured Physical components ACE
vmagahist = histogram(avmag,binsize=20,location=vmagabins,min=-1000,max=1000)
format_hist,vmagahist,vmagabins,vxsize
tahist = histogram(at,binsize=thsize,location=tabins,min=-1000,max=1000)
format_hist,tahist,tabins,thsize
denahist = histogram(aden,binsize=desize,location=denabins,min=-1000,max=1000)
format_hist,denahist,denabins,desize

;plot measured velocity components
;plot measured velocity components
yrange2 = [0.,25.]
plot,vmagdbins,vmagdhist,psym=10, /NOERASE,$
    xtitle='Vmag [km/s]',ytitle='Occurence [%]',position=plot4,$
    xrange=[200.,800.],yrange=yrange2,ystyle=1
    oplot,vmagabins,vmagahist,psym=10,color=200

plot,tdbins,tdhist,psym=10, /NOERASE, $
    xtitle='Th. Speed [km/s]',ytitle='Occurence [%]',position=plot5,$
    xrange=[0.,100.],yrange=yrange2,ystyle=1
    oplot,tabins,tahist,psym=10,color=200

plot,dendbins,dendhist,psym=10, /NOERASE, $
    xtitle='Den. [cm^-3]',ytitle='Occurence [%]',position=plot6,$
    xrange=[0,20.],yrange=yrange2,ystyle=1
    oplot,denabins,denahist,psym=10,color=200

device,/close
;write_png,'../compare_plots/compare_phys_grid.png',tvrd(/true)


;----------------------------------------------------------------------------
;Find  dT/T as a function of V,den,temp to try and remove Th. speed tail
;----------------------------------------------------------------------------
p1 = [.1,.1,.33,.8]
p2 = [.38,.1,.61,.8]
p3 = [.66,.1,.9,.8]
;calc running medians
running_med,wvmag,ddwt,medxdvmag,medydvmag,sig=sigmdvmag,bins=50.
running_med,wt,ddwt,medxdt,medydt,sig=sigmdt,bins=10.
running_med,wden,ddwt,medxdden,medydden,sig=sigmdden,bins=1.5

;Create 2D histograms
dvmag_hist = HIST_2D(wvmag,ddwt,bin1=50.,bin2=1,min1=200.,max1=1000.,min2=-100.,max2=150.)
dt_hist = HIST_2D(wt,ddwt,bin1=10.,bin2=1,min1=0.,max1=200.,min2=-100.,max2=150.)
dden_hist = HIST_2D(wden,ddwt,bin1=3,bin2=1,min1=0.,max1=50.,min2=-100.,max2=150.)

;take the log number
dvmag_hist = alog10(dvmag_hist)
dt_hist = alog10(dt_hist)
dden_hist = alog10(dden_hist)

;replace bad with min.
dvmag_hist[where(finite(dvmag_hist) eq 0)] = min(dvmag_hist,/nan)
dt_hist[where(finite(dt_hist) eq 0)] = min(dt_hist,/nan)
dden_hist[where(finite(dden_hist) eq 0)] = min(dden_hist,/nan)

minv = min(dt_hist)
maxv = max(dt_hist)

;Plot 2d histograms with medians overplotted
cgimage,dvmag_hist ,position=p1,xrange=[200,1000],yrange=[-100.,150.],maxval=maxv,minval=minv
plot,wvmag,ddwt,psym=6,/NODATA,/NOERASE, $
    xtitle="Wind Speed [km/s]",ytitle=delta+"Th. Speed (DSCVOR-WIND)/WIND [%]", position=p1, $
    xrange=[200,1000],yrange=[-100.,150.],xstyle=1,ystyle=1
    oplot,medxdvmag,medydvmag,color=255,psym=7,symsize=1.3
    oplot,medxdvmag,medydvmag,color=0,psym=7
    error_bars,medxdvmag,medydvmag,sigmdvmag

cgimage,dt_hist ,position=p2,xrange=[0.,200.],yrange=[-100.,150.],/NOERASE,maxval=maxv,minval=minv
plot,wt,ddwt,psym=6, /NOERASE, /NODATA, $
    xtitle="Th. Speed [km/s]", position=p2, $
    xrange=[0.,200.],yrange=[-100.,150.],xstyle=1,ystyle=1
    oplot,medxdt,medydt,color=255,psym=7,symsize=1.3
    oplot,medxdt,medydt,color=0,psym=7
    error_bars,medxdt,medydt,sigmdt

cgimage,dden_hist ,position=p3,xrange=[0.,50.],yrange=[-100.,150.],/NOERASE,maxval=maxv,minval=minv
plot,wden,ddwt,psym=6, /NOERASE, /NODATA, $
    xtitle="Density [cm^-3]", position=p3, $
    xrange=[0.,50.],yrange=[-100.,150.],xstyle=1,ystyle=1
    oplot,medxdden,medydden,color=255,psym=7,symsize=1.3
    oplot,medxdden,medydden,color=0,psym=7
    error_bars,medxdden,medydden,sigmdden

cgColorbar,position=[.93,.1,.99,.9],title='log(Number)',textthick=3,charsize=2.0,range=[minv,maxv]


;write_png,'../compare_plots/reg_th_speed_on_X.png',tvrd(/true)


;----------------------------------------------------------------------------
;Find  dT as a function of V,den,temp to try and remove Th. speed tail
;----------------------------------------------------------------------------
;change back to total value
ddwt = dt*ddwt/100.
;calc running medians
running_med,wvmag,ddwt,medxdvmag,medydvmag,sig=sigmdvmag,bins=20.
running_med,wt,ddwt,medxdt,medydt,sig=sigmdt,bins=5.
running_med,wden,ddwt,medxdden,medydden,sig=sigmdden,bins=1.5

;Create 2D histograms
dvmag_hist = HIST_2D(wvmag,ddwt,bin1=40.,bin2=1,min1=200.,max1=1000.,min2=-100.,max2=150.)
dt_hist = HIST_2D(wt,ddwt,bin1=10.,bin2=1,min1=0.,max1=200.,min2=-100.,max2=150.)
dden_hist = HIST_2D(wden,ddwt,bin1=3,bin2=1,min1=0.,max1=50.,min2=-100.,max2=150.)

;take the log number
dvmag_hist = alog10(dvmag_hist)
dt_hist = alog10(dt_hist)
dden_hist = alog10(dden_hist)

;replace bad with min.
dvmag_hist[where(finite(dvmag_hist) eq 0)] = min(dvmag_hist,/nan)
dt_hist[where(finite(dt_hist) eq 0)] = min(dt_hist,/nan)
dden_hist[where(finite(dden_hist) eq 0)] = min(dden_hist,/nan)

minv = min(dt_hist)
maxv = max(dt_hist)

;Plot 2d histograms with medians overplotted
cgimage,dvmag_hist ,position=p1,xrange=[200,1000],yrange=[-100.,150.],maxval=maxv,minval=minv
plot,wvmag,ddwt,psym=6,/NODATA,/NOERASE, $
    xtitle="Wind Speed [km/s]",ytitle=delta+"Th. Speed (DSCVOR-WIND) [km/s]", position=p1, $
    xrange=[200,1000],yrange=[-100.,150.],xstyle=1,ystyle=1
    oplot,medxdvmag,medydvmag,color=255,psym=7,symsize=1.3
    oplot,medxdvmag,medydvmag,color=0,psym=7
    error_bars,medxdvmag,medydvmag,sigmdvmag

cgimage,dt_hist ,position=p2,xrange=[0.,200.],yrange=[-100.,150.],/NOERASE,maxval=maxv,minval=minv
plot,wt,ddwt,psym=6, /NOERASE, /NODATA, $
    xtitle="Th. Speed [km/s]", position=p2, $
    xrange=[0.,200.],yrange=[-100.,150.],xstyle=1,ystyle=1
    oplot,medxdt,medydt,color=255,psym=7,symsize=1.3
    oplot,medxdt,medydt,color=0,psym=7
    error_bars,medxdt,medydt,sigmdt

cgimage,dden_hist ,position=p3,xrange=[0.,50.],yrange=[-100.,150.],/NOERASE,maxval=maxv,minval=minv
plot,wden,ddwt,psym=6, /NOERASE, /NODATA, $
    xtitle="Density [cm^-3]", position=p3, $
    xrange=[0.,50.],yrange=[-100.,150.],xstyle=1,ystyle=1
    oplot,medxdden,medydden,color=255,psym=7,symsize=1.3
    oplot,medxdden,medydden,color=0,psym=7
    error_bars,medxdden,medydden,sigmdden

cgColorbar,position=[.93,.1,.99,.9],title='log(Number)',textthick=3,charsize=2.0,range=[minv,maxv]


;write_png,'../compare_plots/reg_abs_th_speed_on_X.png',tvrd(/true)




end
