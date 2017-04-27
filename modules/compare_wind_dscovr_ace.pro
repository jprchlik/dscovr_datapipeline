;----------------------------------------------------------------------------------------
;
;
;USAGE 
;compare_wind_dscovr_ace,syear,sdoy,edoy=edoy,aceoff=aceoff,eyear=eyear
;
;
;----------------------------------------------------------------------------------------

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

;do the same for ace
        load_plasma,ayear[i],tsdoy,tedoy,/ACE,doy=adoy,bx=abx,by=aby,bz=abz,bmag=abmag,vmag=avmag,vy=avy,vx=avx,vz=avz,t=at,den=aden,x=ax,y=ay,z=az
        load_plasma,ayear[i],tsdoy,tedoy,/WIND,doy=idoy,bx=ibx,by=iby,bz=ibz,bmag=ibmag,vmag=ivmag,vy=ivy,vx=ivx,vz=ivz,t=it,den=iden,x=ix,y=iy,z=iz
    endif else begin
         tsdoy = 1
         load_plasma,ryear[i],tsdoy,tedoy,/DSC,doy=tddoy,bx=tdbx,by=tdby,bz=tdbz,bmag=tdbmag,vmag=tdvmag,vy=tdvy,vx=tdvx,vz=tdvz,t=tdt,den=tdden,x=tdx,y=tdy,z=tdz
         load_plasma,ryear[i],tsdoy,tedoy,/WIND,doy=twdoy,bx=twbx,by=twby,bz=twbz,bmag=twbmag,vmag=twvmag,vy=twvy,vx=twvx,vz=twvz,t=twt,den=twden,x=twx,y=twy,z=twz

;do the same for ace
         load_plasma,ayear[i],tsdoy,tedoy,/ACE,doy=tadoy,bx=tabx,by=taby,bz=tabz,bmag=tabmag,vmag=tavmag,vy=tavy,vx=tavx,vz=tavz,t=tat,den=taden,x=tax,y=tay,z=taz
         load_plasma,ayear[i],tsdoy,tedoy,/WIND,doy=tidoy,bx=tibx,by=tiby,bz=tibz,bmag=tibmag,vmag=tivmag,vy=tivy,vx=tivx,vz=tivz,t=tit,den=tiden,x=tix,y=tiy,z=tiz

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


;Differences in parameters
ddwbx = dbx-wbx
ddwby = dby-wby
ddwbz = dbz-wbz
ddwvx = dvx-wvx
ddwvy = dvy-wvy
ddwvz = dvz-wvz
ddwx = dx-wx
ddwy = dy-wy
ddwz = dz-wz
ddwbmag = dbmag-wbmag
ddwvmag = dvmag-wvmag
ddwt = dt-wt
ddwden = dden-wden

;Differences in parameters
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

plot,vydwbins,vydwhist,psym=10,/NOERASE, $
    xtitle=delta+'Vy (X-WIND) [km/s]',ytitle='Occurence [%]',position=[.41,.6,.65,.95], $
    xrange=vxrange,yrange=vyrange
    oplot,vyaibins,100.*vyaihist/float(total(vyaihist)),psym=10,color=200

plot,vzdwbins,vzdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Vz (X-WIND) [km/s]',ytitle='Occurence [%]',position=[.74,.6,.98,.95], $
    xrange=vxrange,yrange=vyrange
    oplot,vzaibins,vzaihist,psym=10,color=200


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
plot,tdwbins,tdwhist,psym=10,/NOERASE, $
    xtitle=delta+'Th. Speed (X-WIND) [km/s]',ytitle='Occurence [%]',position=[.41,.10,.65,.45], $
    xrange=[-50,50],yrange=vyrange
    oplot,taibins,taihist,psym=10,color=200
plot,dendwbins,dendwhist,psym=10,/NOERASE, $
    xtitle=delta+'Den. (X-WIND) [cm^-3]',ytitle='Occurence [%]',position=[.74,.10,.98,.45], $
    xrange=[-5,5],yrange=vyrange
    oplot,denaibins,denaihist,psym=10,color=200



write_png,'../compare_plots/compare_grid.png',tvrd(/true)


end
