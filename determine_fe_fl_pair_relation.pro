;Determine the relationship between Flares and filament eruptions spatially/temporally separated. 
;Is there a sympathetic relationship?


;----------------------------------------------------------------------------->
;----------------------------------------------------------------------------->
;Get difference in time, dist between subsequent FE/FL in different locations
;1 is the reference, 2 is the events being matched

function determine_fe_fl_pair_relation_link_events_pre_post, tims1, tims2, lon1, lon2, lat1, lat2

minpairdist=20. ;minimum degrees separation between events

;----------------------------------------------------------------------------->
;LOOP OVER REFERENCE EVENTS

n1=n_elements(tims1)
n2=n_elements(tims2)

;Initialise variables
;Preceding flare
tdiff=0.
distdiff=0.
tdifft1=0.
tdifft2=0.

;Following flare
ftdiff=0.
fdistdiff=0.
ftdifft1=0.
ftdifft2=0.

for i=0,n1-1 do begin

;Ref Eruption times
	thist=tims1[i]

;Ref positions
	thislon=theta_shift(tim2carr(thist, offset=lon1[i]))
	thislat=lat1[i]

;----------------------------------------------------------------------------->
;PRECEDING FLARE TO FILAMENT ERUPTION

;NEED TO TAKE FLARES OUTSIDE OF FILAMENT ERUPTION to compute sympathetic candidate waiting times.
;Also limit to >C-class?

	fefldist=gc_dist( rebin([thislon,thislat],2,n2) , transpose([[theta_shift( tim2carr(tims2,offset=lon2) )], [lat2] ]) )

;Find the preceding flare event to each FE but the distance between the two has to be >20deg
	wadjpre=min(where(tims2 gt thist and fefldist gt minpairdist))


;Determine the distances from the big flare to the other flares near in time
	if wadjpre[0] ne -1 then begin

;DEBUG!!!!!!!!!!!!!!!!
;if abs(thist-tims2[wadjpre]) le 3.*3600. then begin & print,anytim(thist,/vms),' ',anytim(tims2[wadjpre],/vms) & stop & endif


		tdiff=[tdiff,thist-tims2[wadjpre]]
		distdiff=[distdiff,gc_dist([thislon,thislat], $
				[theta_shift( tim2carr(tims2[wadjpre],offset=lon2[wadjpre])),lat2[wadjpre] ] $
			)]
		tdifft1=[tdifft1,thist]
		tdifft2=[tdifft2,tims2[wadjpre]]
	endif

;----------------------------------------------------------------------------->
;FOLLOWING FLARE TO FILAMENT ERUPTION

;Find the following flare event to each FE
	wadjpost=max(where(tims2 lt thist and fefldist gt minpairdist))

;Determine the distances from the big flare to the other flares near in time
	if wadjpost[0] ne -1 then begin

		ftdiff=[ftdiff,thist-tims2[wadjpost]]
		fdistdiff=[fdistdiff,gc_dist([thislon,thislat], $
				[theta_shift( tim2carr(tims2[wadjpost],offset=lon2[wadjpost])),lat2[wadjpost] ] $
			)]
		ftdifft1=[ftdifft1,thist]
		ftdifft2=[ftdifft2,tims2[wadjpost]]
	endif

endfor

;preceding flare variables
tdiff=tdiff[1:*]
distdiff=distdiff[1:*]

;following flare variables
ftdiff=ftdiff[1:*]
fdistdiff=fdistdiff[1:*]

tdifft1=tdifft1[1:*]
ftdifft1=ftdifft1[1:*]
tdifft2=tdifft2[1:*]
ftdifft2=ftdifft2[1:*]

outstruct={tdiff:tdiff,t1:tdifft1,t2:tdifft2,distdiff:distdiff,ftdiff:ftdiff,ft1:ftdifft1,ft2:ftdifft2,fdistdiff:fdistdiff}

return,outstruct

end


;----------------------------------------------------------------------------->
;----------------------------------------------------------------------------->
;Link the filament list to the filament eruption list. Get better positions for FEs.

function determine_fe_fl_pair_relation_combine_fi_fe, timfe, timfi, lonfe, lonfi, latfe, latfi, strfe, strfi


;find filament obs just before each filament eruption.
;match in position

;Determine if there was a dissappearance 
;Can the filament found be matched to another filament after the eruption?

;output matched subscripts of FE/FI lists


hcxfe=strmid(strfe.HRC_COORD,6,17)
hcyfe=strmid(strfe.HRC_COORD,24,17)

hcxfi=strmid(strfi.HRC_COORD,6,17)
hcyfi=strmid(strfi.HRC_COORD,24,17)







end


;----------------------------------------------------------------------------->
;----------------------------------------------------------------------------->


pro determine_fe_fl_pair_relation,res1tore=res1tore,res2tore=res2tore

root='~/science/projects/stereo_sympathetic_flaring/'
datapath=root+'data/'
plotpath=root+'plots/determine_fe-fl_pair_relation/png/'
epspath=root+'plots/determine_fe-fl_pair_relation/eps/'
;epspath=root+'eps/'

;Read in the event lists------------------------------------------------------>

if keyword_set(res1tore) then goto,skipres1tore

;All HEK events
restore,datapath+'compile_hek_events-fl_fi_fe_sp_fa.sav',/ver

;HEK flares, separated by wavelength
restore,datapath+'compile_hek_events-flaredetect.sav',/ver

;STEREO 360 Sun list
;restore,datapath+'combine_euvi_events.sav',/ver
restore,datapath+'determine_usual_flare_wait_times-flarelists.sav',/ver

;Last events flares
restore,datapath+'get_lastevents_list.sav',/ver

;Rhessi flares
restore,datapath+'get_rhessi_flares.sav',/ver

restore,datapath+'combine_goesxray_halpha_flare_lists-res1.sav',/ver

;Sort in time and filter the lists for time greater than 2010-04-------------->

;tfilt=anytim('1-Apr-2010 00:00')
tran=[anytim('16-jan-2011'),anytim('11-jun-2013')] ;time range to crop latest events have coverage to match FEs

;Last Events---------------->
nle=n_elements(lescat)
timle=anytim(lescat.date_obs)
stimle=sort(timle)
timle=timle[stimle]
lescat=lescat[stimle]

smart_nsew2hg, lescat.HELIO, latle, lonle
;laslatlon=arcmin2hel(lescat.xcen/60.,lescat.ycen/60.,date=lestest.fpeak)

wle=where(timle ge tran[0] and timle le tran[1])

;Make a list of the flares that are GE M5 class
flrintle=goes_class2int(lescat.class)
;Only for flares occuring after TFILT
wc1=where(flrintle[wle] ge 1d-6)
wc5=where(flrintle[wle] ge 5d-6)
wm1=where(flrintle[wle] ge 1d-5)
wm5=where(flrintle[wle] ge 5d-5)

;----------------------------------------------------------------------------->
goto,skipotherinst ;---------------------------------------------------------->

;Rhessi--------------------->
nrh=n_elements(rhessi_flrs)
timrh=anytim(rhessi_flrs.peak_time)
stimrh=sort(timrh)
timrh=timrh[stimrh]
rhessi_flrs=rhessi_flrs[stimrh]

rhlonlat=fltarr(2,nrh)

;Break the rhessi list in to time blocks of 1 month to convert the HC positions to HG
rhtgrid=anytim(timegrid(anytim(min(timrh))-35*24.*3600.,anytim(max(timrh))+35*24.*3600.,/month))
nrhgrid=n_elements(rhtgrid)
offlimbrh=fltarr(nrhgrid)
for i=0l,nrhgrid-2l do begin
	wthisrhset=where(timrh ge rhtgrid[i] and timrh lt rhtgrid[i+1])
	if wthisrhset[0] ne -1 then begin
		rhlonlat[*,wthisrhset]=arcmin2hel((rhessi_flrs[wthisrhset].position)[1,*]/60.,(rhessi_flrs[wthisrhset].position)[0,wthisrhset]/60.,date=anytim(timrh[median(wthisrhset)],/vms),off_limb=thisofflimbrh,/sphere)
		offlimbrh[wthisrhset]=thisofflimbrh
	endif

endfor
lonrh=rhlonlat[1,*]
latrh=rhlonlat[0,*]

;Also filter for whether the flares are on the limb or not
wrh=where(timrh ge tfilt and offlimbrh eq 0)

;STA/STB Flares------------->
nsta=n_elements(EUVSTRA304)
timsta=anytim(EUVSTRA304.peak)
stimsta=sort(timsta)
timsta=timsta[stimsta]
EUVSTRA304=EUVSTRA304[stimsta]

lonsta=EUVSTRA304.HGLON
latsta=EUVSTRA304.HGLAT

wsta=where(timsta ge tfilt)

nstb=n_elements(EUVSTRB304)
timstb=anytim(EUVSTRB304.peak)
stimstb=sort(timstb)
timstb=timstb[stimstb]
EUVSTRB304=EUVSTRB304[stimstb]

lonstb=EUVSTRB304.HGLON
latstb=EUVSTRB304.HGLAT

wstb=where(timstb ge tfilt)

;HEK Flares----------------->
hekflr=[HEKFLR131,HEKFLR171,HEKFLR193,HEKFLR211,HEKFLR304,HEKFLR335,HEKFLR94]

nhek=n_elements(hekflr)
timhek=anytim(hekflr.EVENT_PEAKTIME)
stimhek=sort(timhek)
timhek=timhek[stimhek]
hekflr=hekflr[stimhek]

heklat=hekflr.HGS_Y
heklon=hekflr.HGS_X
hekclon=hekflr.HGC_X

whek=where(timhek ge tfilt)

;Filament Observations------>

;From Automated method of some kind?? 
;H-alpha data
;Includes length of filament: FI_LENGTH
waafdcc=where(HEKFILAMENT.FRM_NAME eq 'AAFDCC') 
HEKFILAMENT=HEKFILAMENT[waafdcc]

nfi=n_elements(HEKFILAMENT)
timfi=anytim(HEKFILAMENT.EVENT_STARTTIME)
stimfi=sort(timfi)
timfi=timfi[stimfi]
HEKFILAMENT=HEKFILAMENT[stimfi]

lonfi=HEKFILAMENT.HGS_X
latfi=HEKFILAMENT.HGS_Y
clonfi=theta_shift(HEKFILAMENT.HGC_X)

wfi=where(timfi ge tfilt)

;----------------------------------------------------------------------------->
skipotherinst: ;---------------------------------------------------------->

;HEK Filament Eruptions----->
nfe=n_elements(HEKFIERUPT)
timfe=anytim(HEKFIERUPT.EVENT_STARTTIME)
timfeavg=(anytim(HEKFIERUPT.EVENT_STARTTIME)+anytim(HEKFIERUPT.EVENT_ENDTIME))/2.
stimfe=sort(timfe)
timfe=timfe[stimfe]
HEKFIERUPT=HEKFIERUPT[stimfe]

ufe=uniq(HEKFIERUPT.EVENT_STARTTIME+'-'+HEKFIERUPT.EVENT_endtime+'-'+HEKFIERUPT.HGS_COORD)

timfe=timfe[ufe]
HEKFIERUPT=HEKFIERUPT[ufe]

lonfe=float(HEKFIERUPT.HGS_X)
latfe=float(HEKFIERUPT.HGS_Y)
clonfe=theta_shift(float(HEKFIERUPT.HGC_X))

;need to loop over months to make a list of HCXY positions
;hcxy=hel2arcmin(latfe,lonfe,date=)

;stop

wfe=where(finite(timfe) eq 1)
wfe2=where(timfe[wfe] gt tran[0] and timfe[wfe] lt tran[1])


save,file=datapath+'determine_fe-fl_pair_relation.res1.sav'
skipres1tore:
if keyword_set(res1tore) then restore,datapath+'determine_fe-fl_pair_relation.res1.sav',/ver



;--Make a context plot of the event occurrences------------------------------->
setplotenv,file=epspath+'determine_superposed_epoch_flaring-filament_positions_hg_ondisk_andvs_t.eps',/ps,xs=15,ys=5
!p.multi=0
setcolors,/sys,/sil
!p.color=0
!p.background=255
!p.charsize=1.4
!p.thick=5
;Plot event positions in lat,lon
plotsym,0,0.5,/fill
plot,findgen(100),/nodat,/xsty,xran=[-90,90],/ysty,yran=[-90,90],xtit='Longitude',ytit='Latitude',ymargin=[2,2],pos=[0.06,0.16,0.31,0.95] ;,/iso
;oplot,HEKFLR131.hgs_x,HEKFLR131.hgs_y,ps=2,color=!blue
oplot,lonle[wle],latle[wle],ps=1,color=0
oplot,lonfe[wfe[wfe2]],latfe[wfe[wfe2]],ps=8,color=!forest
oplot,lonle[wle[wm1]],latle[wle[wm1]],ps=8,color=!red
;Plot times and flare classes
utplot,timle[wle]-tran[0],flrintle[wle],tran[0],/ylog,ps=3,xtit='Time (2011-2013)',ytit='Flare Class',ymargin=[2,2],/noerase,pos=[0.38,0.16,0.95,0.95],/xsty
;vline,timfe[wfe]-min(timle[wle]),color=!forest,/vlog,yran=[-8,-7]
oplot,timle[wle]-tran[0],flrintle[wle],ps=8,color=0
oplot,timle[wle[wm1]]-tran[0],flrintle[wle[wm1]],ps=8,color=!red
fetimhisty=histogram(timfe[wfe[wfe2]],bin=3600.*24.*7,loc=fetimhistx)
oplot,fetimhistx-tran[0],fetimhisty/1d8,ps=10,color=!forest,thick=5
axis,yaxis=1,yran=[1d-8,1d-3]*1d8,/ylog,color=!forest,chars=1.4,ytit='# FE'
;polyfill,[0.16,0.255,0.255,0.16]-0.09,[0.42,0.42,0.34,0.34]+0.5,color=255,/norm
legend,['All Flares','>M1','FEs'],color=[0,!red,!forest],psym=[4,4,4],/top,/left,chars=1.4,/norm,pos=[0.16-0.09,0.42+0.5]
closeplotenv
;window_Capture,file=plotpath+'determine_superposed_epoch_flaring-filament_positions_hg_ondisk_andvs_t'

stop



;Test C1, M1, M5 FL-FE pairs
;Filament is reference: pre=filament before flare; post=filament after flare
if not keyword_set(res2tore) then begin

strprepostfec1=determine_fe_fl_pair_relation_link_events_pre_post(timfe[wfe[wfe2]], timle[wle[wc1]], lonfe[wfe[wfe2]], lonle[wle[wc1]], latfe[wfe[wfe2]], latle[wle[wc1]])

strprepostfec5=determine_fe_fl_pair_relation_link_events_pre_post(timfe[wfe[wfe2]], timle[wle[wc5]], lonfe[wfe[wfe2]], lonle[wle[wc5]], latfe[wfe[wfe2]], latle[wle[wc5]])

strprepostfem1=determine_fe_fl_pair_relation_link_events_pre_post(timfe[wfe[wfe2]], timle[wle[wm1]], lonfe[wfe[wfe2]], lonle[wle[wm1]], latfe[wfe[wfe2]], latle[wle[wm1]])

strprepostfem5=determine_fe_fl_pair_relation_link_events_pre_post(timfe[wfe[wfe2]], timle[wle[wm5]], lonfe[wfe[wfe2]], lonle[wle[wm5]], latfe[wfe[wfe2]], latle[wle[wm5]])

;Flare is reference: pre=Flare before FE; post=Flare after FE
strprepostflc1=determine_fe_fl_pair_relation_link_events_pre_post(timle[wle[wc1]], timfe[wfe[wfe2]], lonle[wle[wc1]], lonfe[wfe[wfe2]], latle[wle[wc1]], latfe[wfe[wfe2]])

strprepostflc5=determine_fe_fl_pair_relation_link_events_pre_post(timle[wle[wc5]], timfe[wfe[wfe2]], lonle[wle[wc5]], lonfe[wfe[wfe2]], latle[wle[wc5]], latfe[wfe[wfe2]])

strprepostflm1=determine_fe_fl_pair_relation_link_events_pre_post(timle[wle[wm1]], timfe[wfe[wfe2]], lonle[wle[wm1]], lonfe[wfe[wfe2]], latle[wle[wm1]], latfe[wfe[wfe2]])

strprepostflm5=determine_fe_fl_pair_relation_link_events_pre_post(timle[wle[wm5]], timfe[wfe[wfe2]], lonle[wle[wm5]], lonfe[wfe[wfe2]], latle[wle[wm5]], latfe[wfe[wfe2]])

save,strprepostfec1,strprepostfec5,strprepostfem1,strprepostfem5,strprepostflc1,strprepostflc5,strprepostflm1,strprepostflm5,file=datapath+'determine_fe-fl_pair_relation_res1.sav'

endif else restore,datapath+'determine_fe-fl_pair_relation_res1.sav',/ver

stop

;Plot the difference between each pair---------------------------------------->

tplotran=[0.,24.] ;limit pairs to +-24
tbinran=[0.,7.*24.]

;{tdiff:tdiff,distdiff:distdiff,ftdiff:ftdiff,fdistdiff,fdistdiff}

;Filament Eruption is the reference------------->

;W.Time distributions
yhistfeprem5=histogram(-strprepostfem5.tdiff/3600.,bin=3.,loc=xhistfeprem5,min=tplotran[0],max=tbinran[1])
yhistfeprem1=histogram(-strprepostfem1.tdiff/3600.,bin=3.,loc=xhistfeprem1,min=tplotran[0],max=tbinran[1])
yhistfeprec5=histogram(-strprepostfec5.tdiff/3600.,bin=3.,loc=xhistfeprec5,min=tplotran[0],max=tbinran[1])
yhistfeprec1=histogram(-strprepostfec1.tdiff/3600.,bin=3.,loc=xhistfeprec1,min=tplotran[0],max=tbinran[1])

;Pair Distance distributions 
ydistfeprem5=histogram(strprepostfem5.distdiff,bin=1.,loc=xdistfeprem5,min=0,max=180)
ydistfeprem1=histogram(strprepostfem1.distdiff,bin=1.,loc=xdistfeprem1,min=0,max=180)
ydistfeprec5=histogram(strprepostfec5.distdiff,bin=1.,loc=xdistfeprec5,min=0,max=180)
ydistfeprec1=histogram(strprepostfec1.distdiff,bin=1.,loc=xdistfeprec1,min=0,max=180)

;W.Time distributions
yhistfepostm5=histogram(strprepostfem5.ftdiff/3600.,bin=3.,loc=xhistfepostm5,min=tplotran[0],max=tbinran[1])
yhistfepostm1=histogram(strprepostfem1.ftdiff/3600.,bin=3.,loc=xhistfepostm1,min=tplotran[0],max=tbinran[1])
yhistfepostc5=histogram(strprepostfec5.ftdiff/3600.,bin=3.,loc=xhistfepostc5,min=tplotran[0],max=tbinran[1])
yhistfepostc1=histogram(strprepostfec1.ftdiff/3600.,bin=3.,loc=xhistfepostc1,min=tplotran[0],max=tbinran[1])

;Pair Distance distributions 
ydistfepostm5=histogram(strprepostfem5.fdistdiff,bin=1.,loc=xdistfepostm5,min=0,max=180)
ydistfepostm1=histogram(strprepostfem1.fdistdiff,bin=1.,loc=xdistfepostm1,min=0,max=180)
ydistfepostc5=histogram(strprepostfec5.fdistdiff,bin=1.,loc=xdistfepostc5,min=0,max=180)
ydistfepostc1=histogram(strprepostfec1.fdistdiff,bin=1.,loc=xdistfepostc1,min=0,max=180)

;Uncertainties of each distribution assuming poisson
uyhistfeprem5=sqrt(yhistfeprem5)
uyhistfepostm5=sqrt(yhistfepostm5)
uyhistfeprem1=sqrt(yhistfeprem1)
uyhistfepostm1=sqrt(yhistfepostm1)
uyhistfeprec5=sqrt(yhistfeprec5)
uyhistfepostc5=sqrt(yhistfepostc5)
uyhistfeprec1=sqrt(yhistfeprec1)
uyhistfepostc1=sqrt(yhistfepostc1)

xhistfepre=xhistfeprem5 & xhistfepost=xhistfepostm5
;save a structure of the histograms
obsstrfe={yhistfeprem5:yhistfeprem5, yhistfepostm5:yhistfepostm5, $
yhistfeprem1:yhistfeprem1, yhistfepostm1:yhistfepostm1, $
yhistfeprec5:yhistfeprec5, yhistfepostc5:yhistfepostc5, $
yhistfeprec1:yhistfeprec1, yhistfepostc1:yhistfepostc1, $
xhistfepre:xhistfepre, xhistfepost:xhistfepost, $
strprepostfem5:strprepostfem5, strprepostfem1:strprepostfem1, strprepostfec5:strprepostfec5, strprepostfec1:strprepostfec1}
save,obsstrfe,file=datapath+'determine_fe-fl_pair_relation_festruct.sav'

stop

!p.thick=5

;Plot >M5
setplotenv,file=epspath+'determine_fe-fl_pair_relation-compare_fl_fe_pre_post_gtm5_byfe.eps',/ps
!p.charsize=1.4
!p.multi=[0,1,2]
plot,xhistfepostm5,yhistfepostm5,ps=10,yran=[0,20],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfeprem5,yhistfeprem5,ps=10,color=180
plot,xhistfepostm5,yhistfepostm5-yhistfeprem5,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-20,20],/ysty,xtit='M5-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistfepostm5, yhistfepostm5-yhistfeprem5, fltarr(n_elements(xhistfepostm5)), (uyhistfepostm5+uyhistfeprem5),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fe-fl_pair_relation-compare_fl_fe_pre_post_gtm5_byfe'

stop

;Plot >M1
setplotenv,file=epspath+'determine_fe-fl_pair_relation-compare_fl_fe_pre_post_gtm1_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostm1,yhistfepostm1,ps=10,yran=[0,70],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfeprem1,yhistfeprem1,ps=10,color=180
plot,xhistfepostm1,yhistfepostm1-yhistfeprem1,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-60,60],/ysty,xtit='M1-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistfepostm1, yhistfepostm1-yhistfeprem1, fltarr(n_elements(xhistfepostm1)), (uyhistfepostm1+uyhistfeprem1),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fe-fl_pair_relation-compare_fl_fe_pre_post_gtm1_byfe'

stop

;Plot >C5
setplotenv,file=epspath+'determine_fe-fl_pair_relation-compare_fl_fe_pre_post_gtc5_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostc5,yhistfepostc5,ps=10,yran=[0,120],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfeprec5,yhistfeprec5,ps=10,color=180
plot,xhistfepostc5,yhistfepostc5-yhistfeprec5,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-60,60],/ysty,xtit='C5-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistfepostc5, yhistfepostc5-yhistfeprec5, fltarr(n_elements(xhistfepostc5)), (uyhistfepostc5+uyhistfeprec5),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fe-fl_pair_relation-compare_fl_fe_pre_post_gtc5_byfe'

stop

;Plot >C1
setplotenv,file=epspath+'determine_fe-fl_pair_relation-compare_fl_fe_pre_post_gtc1_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostc1,yhistfepostc1,ps=10,yran=[0,600],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfeprec1,yhistfeprec1,ps=10,color=180
plot,xhistfepostc1,yhistfepostc1-yhistfeprec1,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-150,150],/ysty,xtit='C1-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistfepostc1, yhistfepostc1-yhistfeprec1, fltarr(n_elements(xhistfepostc1)), (uyhistfepostc1+uyhistfeprec1),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fe-fl_pair_relation-compare_fl_fe_pre_post_gtc1_byfe'

stop

;Plot the y-axis as the number of events in a given bin / total number of event pairs within 5 days
;Plot >M5, Fractional--------------------------------------------------------->
!p.charsize=2
!p.multi=[0,1,2]
!p.background=255
!p.color=0

stop

npostm5=float(n_elements(strprepostfem5.ftdiff))
npostm5hist=float(total(yhistfepostm5))
nprem5=float(n_elements(strprepostfem5.tdiff))
nprem5hist=float(total(yhistfeprem5))
yhistfeavgm5=(yhistfepostm5+yhistfeprem5)/2.

setplotenv,file=epspath+'determine_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtm5_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostm5,yhistfepostm5/npostm5,ps=10,yran=[0,.02],ytit='#/Ntot (Fl. bef.=blk, aft.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=180
oplot,xhistfeprem5,yhistfeprem5/nprem5,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistfepostm5[1:7])/npostm5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfepostm5[8:*])/npostm5),color=0,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistfeprem5[1:7])/nprem5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfeprem5[8:*])/nprem5),color=180,/norm
plot,xhistfepostm5,(yhistfepostm5-yhistfeprem5)/yhistfeavgm5,xran=[0,24],ytit='(Flare bef. - Flare aft.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='M5-FE Wait-time [hrs]'
vline,3,lines=2,color=180
hline,0,thick=2,color=150
oploterr, xhistfepostm5, (yhistfepostm5-yhistfeprem5)/yhistfeavgm5, fltarr(n_elements(xhistfepostm5)), (uyhistfepostm5+uyhistfeprem5)/yhistfeavgm5,ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtm5_byfe'

stop

;Plot >M1, Fractional
npostm1=float(n_elements(strprepostfem1.ftdiff))
npostm1hist=float(total(yhistfepostm1))
nprem1=float(n_elements(strprepostfem1.tdiff))
nprem1hist=float(total(yhistfeprem1))
yhistfeavgm1=(yhistfepostm1+yhistfeprem1)/2.

setplotenv,file=epspath+'determine_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtm1_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostm1,yhistfepostm1/npostm1,ps=10,yran=[0,.05],ytit='# / Ntot',ymarg=[2,1],xran=[0,24],xtit='',xtickname=strarr(10)+' ',ytickname=[' ']
sharpcorners
oploterr,xhistfepostm1,yhistfepostm1/npostm1,uyhistfepostm1/npostm1,ps=4
legend,['Flare Before','Flare After'],lines=[0,0],color=[0,180],/bottom,/right
;vline,3,lines=2,color=180
oplot,xhistfeprem1,yhistfeprem1/nprem1,ps=10,color=180
!p.color=180
oploterr,xhistfeprem1,yhistfeprem1/nprem1,uyhistfeprem1/nprem1,ps=4
!p.color=0
;xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistfepostm1[1:7])/npostm1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfepostm1[8:*])/npostm1),color=0,/norm
;xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistfeprem1[1:7])/nprem1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfeprem1[8:*])/nprem1),color=180,/norm
plot,xhistfepostm1,(yhistfepostm1-yhistfeprem1)/yhistfeavgm1,ytit='(Flare Before - Flare After) / Avg(Before,After)',ps=10,yran=[-2,2],/ysty,xtit='M1-FE Wait-time [hrs]',xran=[0,24],ymarg=[5,-2]
sharpcorners
;vline,3,lines=2,color=180
hline,0,color=150
oploterr, xhistfepostm1, (yhistfepostm1-yhistfeprem1)/yhistfeavgm1, fltarr(n_elements(xhistfepostm1)), (uyhistfepostm1+uyhistfeprem1)/yhistfeavgm1,ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtm1_byfe'


stop

;Plot >C5, Fractional
npostc5=float(n_elements(strprepostfec5.ftdiff))
npostc5hist=float(total(yhistfepostc5))
nprec5=float(n_elements(strprepostfec5.tdiff))
nprec5hist=float(total(yhistfeprec5))
yhistfeavgc5=(yhistfepostc5+yhistfeprec5)/2.

setplotenv,file=epspath+'determine_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtc5_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostc5,yhistfepostc5/npostc5,ps=10,yran=[0,.08],ytit='#/Ntot (Fl. bef.=blk, aft.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=180
oplot,xhistfeprec5,yhistfeprec5/nprec5,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistfepostc5[1:7])/npostc5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfepostc5[8:*])/npostc5),color=0,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistfeprec5[1:7])/nprec5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfeprec5[8:*])/nprec5),color=180,/norm
plot,xhistfepostc5,(yhistfepostc5-yhistfeprec5)/yhistfeavgc5,ytit='(Flare bef. - Flare aft.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='M5-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=180
hline,0,thick=2,color=150
oploterr, xhistfepostc5, (yhistfepostc5-yhistfeprec5)/yhistfeavgc5, fltarr(n_elements(xhistfepostc5)), (uyhistfepostc5+uyhistfeprec5)/yhistfeavgc5,ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtc5_byfe'

stop

;Plot >C1, Fractional
npostc1=float(n_elements(strprepostfec1.ftdiff))
npostc1hist=float(total(yhistfepostc1))
nprec1=float(n_elements(strprepostfec1.tdiff))
nprec1hist=float(total(yhistfeprec1))
yhistfeavgc1=(yhistfepostc1+yhistfeprec1)/2.

setplotenv,file=epspath+'determine_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtc1_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostc1,yhistfepostc1/npostc1,ps=10,yran=[0,.16],ytit='#/Ntot (Fl. bef.=blk, aft.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=180
oplot,xhistfeprec1,yhistfeprec1/nprec1,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistfepostc1[1:7])/npostc1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfepostc1[8:*])/npostc1),color=0,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistfeprec1[1:7])/nprec1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfeprec1[8:*])/nprec1),color=180,/norm
plot,xhistfepostc1,(yhistfepostc1-yhistfeprec1)/yhistfeavgc1,ytit='(Flare bef. - Flare aft.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='M5-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=180
hline,0,thick=2,color=150
oploterr, xhistfepostc1, (yhistfepostc1-yhistfeprec1)/yhistfeavgc1, fltarr(n_elements(xhistfepostc1)), (uyhistfepostc1+uyhistfeprec1)/yhistfeavgc1,ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtc1_byfe'


stop


;Flare is the reference------------------------->

;W.Time distributions
yhistflprem5=histogram(-strprepostflm5.tdiff/3600.,bin=3.,loc=xhistflprem5,min=tplotran[0],max=tbinran[1])
yhistflprem1=histogram(-strprepostflm1.tdiff/3600.,bin=3.,loc=xhistflprem1,min=tplotran[0],max=tbinran[1])
yhistflprec5=histogram(-strprepostflc5.tdiff/3600.,bin=3.,loc=xhistflprec5,min=tplotran[0],max=tbinran[1])
yhistflprec1=histogram(-strprepostflc1.tdiff/3600.,bin=3.,loc=xhistflprec1,min=tplotran[0],max=tbinran[1])

;Pair Distance distributions 
ydistflprem5=histogram(strprepostflm5.distdiff,bin=1.,loc=xdistflprem5,min=0,max=180)
ydistflprem1=histogram(strprepostflm1.distdiff,bin=1.,loc=xdistflprem1,min=0,max=180)
ydistflprec5=histogram(strprepostflc5.distdiff,bin=1.,loc=xdistflprec5,min=0,max=180)
ydistflprec1=histogram(strprepostflc1.distdiff,bin=1.,loc=xdistflprec1,min=0,max=180)

;W.Time distributions
yhistflpostm5=histogram(strprepostflm5.ftdiff/3600.,bin=3.,loc=xhistflpostm5,min=tplotran[0],max=tbinran[1])
yhistflpostm1=histogram(strprepostflm1.ftdiff/3600.,bin=3.,loc=xhistflpostm1,min=tplotran[0],max=tbinran[1])
yhistflpostc5=histogram(strprepostflc5.ftdiff/3600.,bin=3.,loc=xhistflpostc5,min=tplotran[0],max=tbinran[1])
yhistflpostc1=histogram(strprepostflc1.ftdiff/3600.,bin=3.,loc=xhistflpostc1,min=tplotran[0],max=tbinran[1])

;Pair Distance distributions 
ydistflpostm5=histogram(strprepostflm5.fdistdiff,bin=1.,loc=xdistflpostm5,min=0,max=180)
ydistflpostm1=histogram(strprepostflm1.fdistdiff,bin=1.,loc=xdistflpostm1,min=0,max=180)
ydistflpostc5=histogram(strprepostflc5.fdistdiff,bin=1.,loc=xdistflpostc5,min=0,max=180)
ydistflpostc1=histogram(strprepostflc1.fdistdiff,bin=1.,loc=xdistflpostc1,min=0,max=180)

;Uncertainties of each distribution assuming poisson
uyhistflprem5=sqrt(yhistflprem5)
uyhistflpostm5=sqrt(yhistflpostm5)
uyhistflprem1=sqrt(yhistflprem1)
uyhistflpostm1=sqrt(yhistflpostm1)
uyhistflprec5=sqrt(yhistflprec5)
uyhistflpostc5=sqrt(yhistflpostc5)
uyhistflprec1=sqrt(yhistflprec1)
uyhistflpostc1=sqrt(yhistflpostc1)

xhistflpre=xhistflprem5 & xhistflpost=xhistflpostm5
;save a structure of the histograms
obsstrfl={yhistflprem5:yhistflprem5, yhistflpostm5:yhistflpostm5, $
yhistflprem1:yhistflprem1, yhistflpostm1:yhistflpostm1, $
yhistflprec5:yhistflprec5, yhistflpostc5:yhistflpostc5, $
yhistflprec1:yhistflprec1, yhistflpostc1:yhistflpostc1, $
xhistflpre:xhistflpre, xhistflpost:xhistflpost, $
strprepostflm5:strprepostflm5, strprepostflm1:strprepostflm1, strprepostflc5:strprepostflc5, strprepostflc1:strprepostflc1}
save,obsstrfl,file=datapath+'determine_fe-fl_pair_relation_flstruct.sav'

;Plot >M5

!p.charsize=2
!p.multi=[0,1,2]

setplotenv,file=epspath+'determine_fl-fe_pair_relation-compare_fl_fe_pre_post_gtm5_byfe.eps',/ps
!p.charsize=1.4
!p.multi=[0,1,2]
plot,xhistflprem5,yhistflprem5,ps=10,yran=[0,20],xran=[0,24],ytit='# (FE. aft.=blk, bef.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistflpostm5,yhistflpostm5,ps=10,color=180
plot,xhistflpostm5,yhistflprem5-yhistflpostm5,xran=[0,24],ytit='(FE after Flare) - (FE before Flare)',ps=10,yran=[-20,20],/ysty,xtit='M5-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprem5, yhistflprem5-yhistflpostm5, fltarr(n_elements(xhistflprem5)), (uyhistflprem5+uyhistflpostm5),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fl-fe_pair_relation-compare_fl_fe_pre_post_gtm5_byfe'

stop

;Plot >M1

setplotenv,file=epspath+'determine_fl-fe_pair_relation-compare_fl_fe_pre_post_gtm1_byfe.eps',/ps,thick=5,xthick=7,ythick=7
!p.charsize=1.4
plot,xhistflprem1,yhistflprem1/nprem1,ps=10,yran=[0,.1],ytit='# / Ntot',ymarg=[2,1],xran=[0,24],xtit='',xtickname=strarr(10)+' ',ytickname=[' ']
sharpcorners
oploterr,xhistflprem1,yhistflprem1/npostm1,uyhistflprem1/nprem1,ps=4
legend,['FE Before','FE After'],lines=[0,0],color=[0,180],/bottom,/right
;vline,3,lines=2,color=150
oplot,xhistflpostm1,yhistflpostm1/npostm1,ps=10,color=180
!p.color=180
oploterr,xhistflpostm1,yhistflpostm1/npostm1,uyhistflpostm1/npostm1,ps=4
!p.color=0
;xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistflpostm1[1:7])/npostm1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflpostm1[8:*])/npostm1),color=180,/norm
;xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistflprem1[1:7])/nprem1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflprem1[8:*])/nprem1),color=0,/norm
plot,xhistflpostm1,(yhistflprem1-yhistflpostm1),ytit='(FE After - FE Before) / Avg(Before,After)',ps=10,yran=[-2,2],/ysty,xtit='M1-FE Wait-time [hrs]',xran=[0,24],ymarg=[5,-2]
sharpcorners
;vline,3,lines=2,color=150
hline,0,color=150
oploterr, xhistflprem1, (yhistflprem1-yhistflpostm1), fltarr(n_elements(xhistflprem1)), (uyhistflprem1+uyhistflpostm1),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fl-fe_pair_relation-compare_fl_fe_pre_post_gtm1_byfe'

stop

;Plot >C5

setplotenv,file=epspath+'determine_fl-fe_pair_relation-compare_fl_fe_pre_post_gtc5_byfe.eps',/ps
!p.charsize=1.4
!p.multi=[0,1,2]
plot,xhistflprec5,yhistflprec5,ps=10,xran=[0,24],yran=[0,50],ytit='# (FE. aft.=blk, bef.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistflpostc5,yhistflpostc5,ps=10,color=180
plot,xhistflpostc5,yhistflprec5-yhistflpostc5,xran=[0,24],ytit='(FE after Flare) - (FE before Flare)',ps=10,yran=[-20,20],/ysty,xtit='M5-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprec5, yhistflprec5-yhistflpostc5, fltarr(n_elements(xhistflprec5)), (uyhistflprec5+uyhistflpostc5),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fl-fe_pair_relation-compare_fl_fe_pre_post_gtc5_byfe'

stop

;Plot >C1

setplotenv,file=epspath+'determine_fl-fe_pair_relation-compare_fl_fe_pre_post_gtc1_byfe.eps',/ps
!p.charsize=1.4
!p.multi=[0,1,2]
plot,xhistflprec1,yhistflprec1,ps=10,xran=[0,24],yran=[0,300],ytit='# (FE. aft.=blk, bef.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistflpostc1,yhistflpostc1,ps=10,color=180
plot,xhistflpostc1,yhistflprec1-yhistflpostc1,xran=[0,24],ytit='(FE after Flare) - (FE before Flare)',ps=10,yran=[-50,50],/ysty,xtit='M1-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprec1, yhistflprec1-yhistflpostc1, fltarr(n_elements(xhistflprec1)), (uyhistflprec1+uyhistflpostc1),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fl-fe_pair_relation-compare_fl_fe_pre_post_gtc1_byfe'


stop

;Plot >M5; Fractional--------------------------------------------------------->

npostm5=float(n_elements(strprepostflm5.ftdiff))
npostm5hist=float(total(yhistflpostm5))
nprem5=float(n_elements(strprepostflm5.tdiff))
nprem5hist=float(total(yhistflprem5))
yhistflavgm5=(yhistflpostm5+yhistflprem5)/2.

setplotenv,file=epspath+'determine_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtm5_byfe.eps',/ps
!p.charsize=1.4
plot,xhistflprem5,yhistflprem5/nprem5,ps=10,yran=[0,.2],ytit='#/Ntot (FE. aft.=blk, bef.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=150
oplot,xhistflpostm5,yhistflpostm5/npostm5,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistflpostm5[1:7])/npostm5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflpostm5[8:*])/npostm5),color=180,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistflprem5[1:7])/nprem5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflprem5[8:*])/nprem5),color=0,/norm
plot,xhistflpostm5,(yhistflprem5-yhistflpostm5)/yhistflavgm5,ytit='(FE aft. - FE bef.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='M5-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprem5, (yhistflprem5-yhistflpostm5)/yhistflavgm5, fltarr(n_elements(xhistflprem5)), (uyhistflprem5+uyhistflpostm5)/yhistflavgm5,ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtm5_byfe'

stop

;Plot >M1; Fractional

npostm1=float(n_elements(strprepostflm1.ftdiff))
npostm1hist=float(total(yhistflpostm1))
nprem1=float(n_elements(strprepostflm1.tdiff))
nprem1hist=float(total(yhistflprem1))
yhistflavgm1=(yhistflpostm1+yhistflprem1)/2.

setplotenv,file=epspath+'determine_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtm1_byfe.eps',/ps
!p.charsize=1.4
plot,xhistflprem1,yhistflprem1/nprem1,ps=10,yran=[0,.1],ytit='# / Ntot',ymarg=[2,1],xran=[0,24],xtit='',xtickname=strarr(10)+' ',ytickname=[' ']
sharpcorners
oploterr,xhistflprem1,yhistflprem1/nprem1,uyhistfeprem1/nprem1,ps=4
legend,['FE After','FE Before'],lines=[0,0],color=[0,180],/bottom,/right
;vline,3,lines=2,color=150
oplot,xhistflpostm1,yhistflpostm1/npostm1,ps=10,color=180
!p.color=180
oploterr,xhistflpostm1,yhistflpostm1/npostm1,uyhistflpostm1/npostm1,ps=4
!p.color=0
;xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistflpostm1[1:7])/npostm1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflpostm1[8:*])/npostm1),color=180,/norm
;xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistflprem1[1:7])/nprem1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflprem1[8:*])/nprem1),color=0,/norm
plot,xhistflpostm1,(yhistflprem1-yhistflpostm1)/yhistflavgm1,ytit='(FE after - FE before)/Avg(Before,After)',ps=10,yran=[-2,2],/ysty,xtit='M1-FE Wait-time [hrs]',xran=[0,24],ymarg=[5,-2]
sharpcorners
;vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprem1, (yhistflprem1-yhistflpostm1)/yhistflavgm1, fltarr(n_elements(xhistflprem1)), (uyhistflprem1+uyhistflpostm1)/yhistflavgm1,ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtm1_byfe'

stop

;Plot >C5; Fractional

npostc5=float(n_elements(strprepostflc5.ftdiff))
npostc5hist=float(total(yhistflpostc5))
nprec5=float(n_elements(strprepostflc5.tdiff))
nprec5hist=float(total(yhistflprec5))
yhistflavgc5=(yhistflpostc5+yhistflprec5)/2.

setplotenv,file=epspath+'determine_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtc5_byfe.eps',/ps
!p.charsize=1.4
plot,xhistflprec5,yhistflprec5/nprec5,ps=10,yran=[0,.1],ytit='#/Ntot (FE. aft.=blk, bef.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=150
oplot,xhistflpostc5,yhistflpostc5/npostc5,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistflpostc5[1:7])/npostc5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflpostc5[8:*])/npostc5),color=180,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistflprec5[1:7])/nprec5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflprec5[8:*])/nprec5),color=0,/norm
plot,xhistflpostc5,(yhistflprec5-yhistflpostc5)/yhistflavgc5,ytit='(FE aft. - FE bef.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='C5-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprec5, (yhistflprec5-yhistflpostc5)/yhistflavgc5, fltarr(n_elements(xhistflprec5)), (uyhistflprec5+uyhistflpostc5)/yhistflavgc5,ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtc5_byfe'

stop

;Plot >C1; Fractional

npostc1=float(n_elements(strprepostflc1.ftdiff))
npostc1hist=float(total(yhistflpostc1))
nprec1=float(n_elements(strprepostflc1.tdiff))
nprec1hist=float(total(yhistflprec1))
yhistflavgc1=(yhistflpostc1+yhistflprec1)/2.

setplotenv,file=epspath+'determine_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtc1_byfe.eps',/ps
!p.charsize=1.4
plot,xhistflprec1,yhistflprec1/nprec1,ps=10,yran=[0,.1],ytit='#/Ntot (FE. aft.=blk, bef.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=150
oplot,xhistflpostc1,yhistflpostc1/npostc1,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistflpostc1[1:7])/npostc1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflpostc1[8:*])/npostc1),color=180,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistflprec1[1:7])/nprec1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflprec1[8:*])/nprec1),color=0,/norm
plot,xhistflpostc1,(yhistflprec1-yhistflpostc1)/yhistflavgc1,ytit='(FE aft. - FE bef.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='C5-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprec1, (yhistflprec1-yhistflpostc1)/yhistflavgc1, fltarr(n_elements(xhistflprec1)), (uyhistflprec1+uyhistflpostc1)/yhistflavgc1,ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'determine_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtc1_byfe'

;!!!FIX Karels corrections


;!!!Run simulations to test our hypothesis
;-random/random and random/clumping
;-Pg 46 of lab notebook

stop

end
