;Determine the relationship between Flares and filament eruptions spatially/temporally separated. 
;Is there a sympathetic relationship?


;----------------------------------------------------------------------------->
;----------------------------------------------------------------------------->
;Get difference in time, dist between subsequent FE/FL in different locations
;1 is the reference, 2 is the events being matched

function simulate_fe_fl_pair_relation_link_events_pre_post, tims1, tims2, lon1, lon2, lat1, lat2, noposition=noposition

if not keyword_set(noposition) then noposition=0 else noposition=1

minpairdist=20. ;minimum degrees separation between events

;----------------------------------------------------------------------------->
;LOOP OVER REFERENCE EVENTS

n1=n_elements(tims1)
n2=n_elements(tims2)

;Initialise variables
;Preceding flare
tdiff=0.
distdiff=0.

;Following flare
ftdiff=0.
fdistdiff=0.

for i=0,n1-1 do begin

;Ref Eruption times
	thist=tims1[i]

;Ref positions
	if noposition then begin
		thislon=0.
		thislat=0.
	endif else begin
		thislon=theta_shift(tim2carr(thist, offset=lon1[i]))
		thislat=lat1[i]
	endelse

;----------------------------------------------------------------------------->
;PRECEDING FLARE TO FILAMENT ERUPTION

;NEED TO TAKE FLARES OUTSIDE OF FILAMENT ERUPTION to compute sympathetic candidate waiting times.
;Also limit to >C-class?
	if not noposition then $
		fefldist=gc_dist( rebin([thislon,thislat],2,n2) , transpose([[theta_shift( tim2carr(tims2,offset=lon2) )], [lat2] ]) )

;Find the preceding flare event to each FE but the distance between the two has to be >20deg
	if not noposition then $	
		wadjpre=min(where(tims2 gt thist and fefldist gt minpairdist)) $
	else wadjpre=min(where(tims2 gt thist))


;Determine the distances from the big flare to the other flares near in time
	if wadjpre[0] ne -1 then begin

		tdiff=[tdiff,thist-tims2[wadjpre]]
		if not noposition then $
			distdiff=[distdiff,gc_dist([thislon,thislat], $
				[theta_shift( tim2carr(tims2[wadjpre],offset=lon2[wadjpre])),lat2[wadjpre] ] $
			)]
	endif

;----------------------------------------------------------------------------->
;FOLLOWING FLARE TO FILAMENT ERUPTION

;Find the following flare event to each FE
	if not noposition then $
		wadjpost=max(where(tims2 lt thist and fefldist gt minpairdist)) $
	else wadjpost=max(where(tims2 lt thist))
		

;Determine the distances from the big flare to the other flares near in time
	if wadjpost[0] ne -1 then begin

		ftdiff=[ftdiff,thist-tims2[wadjpost]]
		if not noposition then $
			fdistdiff=[fdistdiff,gc_dist([thislon,thislat], $
				[theta_shift( tim2carr(tims2[wadjpost],offset=lon2[wadjpost])),lat2[wadjpost] ] $
			)]
	endif

endfor

;preceding flare variables
tdiff=tdiff[1:*]
if not noposition then distdiff=distdiff[1:*]

;following flare variables
ftdiff=ftdiff[1:*]
if not noposition then fdistdiff=fdistdiff[1:*]

if not noposition then outstruct={tdiff:tdiff,distdiff:distdiff,ftdiff:ftdiff,fdistdiff:fdistdiff} $
	else outstruct={tdiff:tdiff,ftdiff:ftdiff}

return,outstruct

end


;----------------------------------------------------------------------------->
;----------------------------------------------------------------------------->
;Link the filament list to the filament eruption list. Get better positions for FEs.

function simulate_fe_fl_pair_relation_combine_fi_fe, timfe, timfi, lonfe, lonfi, latfe, latfi, strfe, strfi


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

;Generate Simulated Time Series of FEs and FLs
function simulate_fe_fl_pair_relation_tseries, obst

;Find numbers of events
nc1=n_elements(obst)

t1c1=min(obst)
t2c1=max(obst)
tranc1=t2c1-t1c1

;Determine mean flare rate for each flare class, and the FEs

flareratec1=nc1/(t2c1-t1c1)

;Generate waiting times for each, based on the flare rate

flarewtdistc1=exponential(min=0,max=10,/monte,nmont=nc1,outdist=flarewtc1)

;Generate time series based on those waiting times

flaretc1=flare_wait_time(flarewtc1,/inverse)

;Scale to the same as observations

flaretc1=flaretc1-min(flaretc1)
flaretc1=(flaretc1)/max(flaretc1)*tranc1

return,flaretc1

end

;----------------------------------------------------------------------------->
;----------------------------------------------------------------------------->

;Run simulations of the FE / FL WTs to compare to the observed distributions
pro simulate_fe_fl_pair_relation

root='~/science/projects/stereo_sympathetic_flaring/'
datapath=root+'data/'
plotpath=root+'plots/determine_fe-fl_pair_relation/'
epspath=root+'plots/determine_fe-fl_pair_relation/eps/'
;epspath=root+'eps/'

restore,datapath+'determine_fe-fl_pair_relation.res1.sav',/ver

;1. Simulate Flare occurrence using observed flare rates and assuming poisson distribution
;!!! Do I need to care about the position of the flares/FEs? or just do it completely random.
;!!! Depends what percentage of FEs would have been connected to Flares, compared to all possible...



;;Test C1, M1, M5 FL-FE pairs
;;Filament is reference: pre=filament before flare; post=filament after flare
;strprepostfec1=determine_fe_fl_pair_relation_link_events_pre_post(timfe[wfe[wfe2]], timle[wle[wc1]], ;lonfe[wfe[wfe2]], lonle[wle[wc1]], latfe[wfe[wfe2]], latle[wle[wc1]])

;Observed timeseries
tlec1 = timle[wle[wc1]]
tlec5 = timle[wle[wc5]]
tleM1 = timle[wle[wM1]]
tleM5 = timle[wle[wM5]]
tfe = timfe[wfe[wfe2]]

stop

fratec1=n_elements(tlec1)/(max(tlec1)-min(tlec1))*24.*3600.
fratec5=n_elements(tlec5)/(max(tlec5)-min(tlec5))*24.*3600.
fratem1=n_elements(tlem1)/(max(tlem1)-min(tlem1))*24.*3600.
fratem5=n_elements(tlem5)/(max(tlem5)-min(tlem5))*24.*3600.
fratefe=n_elements(tfe)/(max(tfe)-min(tfe))*24.*3600.

simtlec1=simulate_fe_fl_pair_relation_tseries(tlec1)
simtlec5=simulate_fe_fl_pair_relation_tseries(tlec5)
simtlem1=simulate_fe_fl_pair_relation_tseries(tlem1)
simtlem5=simulate_fe_fl_pair_relation_tseries(tlem5)
simtfe=simulate_fe_fl_pair_relation_tseries(tfe)

simfratec1=n_elements(simtlec1)/(max(simtlec1)-min(simtlec1))*24.*3600.
simfratec5=n_elements(simtlec5)/(max(simtlec5)-min(simtlec5))*24.*3600.
simfratem1=n_elements(simtlem1)/(max(simtlem1)-min(simtlem1))*24.*3600.
simfratem5=n_elements(simtlem5)/(max(simtlem5)-min(simtlem5))*24.*3600.
simfratefe=n_elements(simtfe)/(max(simtfe)-min(simtfe))*24.*3600.

plot_hist,tlec1,bin=3600.*24.*7.,/xsty
plot_hist,simtlec1,bin=3600.*24.*7.,/ysty
plot_hist,histogram(tlec1,bin=3600.*24.*7.),bin=5,xran=[1,100]
plot_hist,histogram(simtlec1,bin=3600.*24.*7.),bin=5,xran=[1,100]

print,'rate c1-obs/sim',fratec1,simfratec1
print,'rate c5-obs/sim',fratec5,simfratec5
print,'rate m1-obs/sim',fratem1,simfratem1
print,'rate m5-obs/sim',fratem5,simfratem5
print,'rate fe-obs/sim',fratefe,simfratefe

stop

;compare the simulated and observed waiting time distributions

ytlem5=histogram(flare_wait_time(tlem5)/3600.,bin=3,loc=xtlem5)
uytlem5=sqrt(ytlem5)
ysimtlem5=histogram(flare_wait_time(simtlem5)/3600.,bin=3,loc=xsimtlem5)
uysimtlem5=sqrt(ysimtlem5)
ytfe=histogram(flare_wait_time(tfe)/3600.,bin=3,loc=xtfe)
uytfe=sqrt(ytfe)
ysimtfe=histogram(flare_wait_time(simtfe)/3600.,bin=3,loc=xsimtfe)
uysimtfe=sqrt(ysimtfe)

setplotenv,file=epspath+'simulate_fe-fl_pair_relation-compare_sim_obs_wts.eps',/ps,ys=7,xs=21
!p.multi=[0,2,1]
!p.background=255
!p.color=0
!p.charsize=2
!x.thick=5
!y.thick=5
!p.thick=4
;plot_hist,flare_wait_time(tlem1)/3600.,xran=[0,96],/log,yran=[0.1,1d2],xtit='Observed Waiting-times [hrs]',ytit='# Event Pairs',bin=3
plot,xtlem5,ytlem5,xran=[0,96],/ylog,/xsty,/ysty,yran=[0.5,1d2],xtit='Observed Waiting-times [hrs]',ytit='# Event Pairs',ps=10,xmargin=[7,-2]
sharpcorners,thick=4
oploterr, xtlem5, ytlem5, uytlem5,ps=4,ymarg=2
;vline,[30,51,69,72,75,81],yran=alog10([0.1,2.]),/vlog
legend,['>M5 Flares','Filament Eruptions'],/top,/right,color=[0,150],lines=[0,0],chars=1.4
sharpcorners
;plot_hist,flare_wait_time(tfe)/3600.,/oplot,color=150,bin=2
oplot,xtfe,ytfe,color=150,ps=10
!p.color=150
oploterr, xtfe, ytfe, uytfe,ps=4,ymarg=2,color=150
vline,90,yran=alog10([0.1,2.]),/vlog
!p.color=0
;plot_hist,flare_wait_time(simtlem1)/3600.,xran=[0,96],/log,yran=[0.1,1d2],xtit='Simulated Waiting-times [hrs]',ytit='# Event Pairs',bin=3
plot,xsimtlem5,ysimtlem5,xran=[0,96],/ylog,/xsty,/ysty,yran=[0.5,1d2],xtit='Simulated Waiting-times [hrs]',ytit='',ytickname=strarr(10)+' ',ps=10,xmargin=[2,3]
sharpcorners,thick=4
oploterr, xsimtlem5, ysimtlem5, uysimtlem5,ps=4,ymarg=2
;vline,[42,66,87,90],yran=alog10([0.1,2.]),/vlog
;plot_hist,flare_wait_time(simtfe)/3600.,/oplot,color=150,bin=2
oplot,xsimtfe,ysimtfe,color=150,ps=10
!p.color=150
oploterr, xsimtfe, ysimtfe, uysimtfe,ps=4,ymarg=2,color=150
!p.color=0
closeplotenv

stop

;Test C1, M1, M5 FL-FE pairs
;Filament is reference: pre=filament before flare; post=filament after flare
strprepostfec1=simulate_fe_fl_pair_relation_link_events_pre_post(simtfe, simtlec1,/nopos)
strprepostfec5=simulate_fe_fl_pair_relation_link_events_pre_post(simtfe, simtlec5,/nopos)
strprepostfem1=simulate_fe_fl_pair_relation_link_events_pre_post(simtfe, simtlem1,/nopos)
strprepostfem5=simulate_fe_fl_pair_relation_link_events_pre_post(simtfe, simtlem5,/nopos)

strprepostflc1=simulate_fe_fl_pair_relation_link_events_pre_post(simtlec1,simtfe,/nopos)
strprepostflc5=simulate_fe_fl_pair_relation_link_events_pre_post(simtlec5,simtfe, /nopos)
strprepostflm1=simulate_fe_fl_pair_relation_link_events_pre_post(simtlem1,simtfe, /nopos)
strprepostflm5=simulate_fe_fl_pair_relation_link_events_pre_post(simtlem5,simtfe, /nopos)

;Plot the difference between each pair---------------------------------------->

tplotran=[0.,24.] ;limit pairs to +-24
tbinran=[0.,7.*24.]

;Filament Eruption is the reference------------->

;W.Time distributions
yhistfeprem5=histogram(-strprepostfem5.tdiff/3600.,bin=3.,loc=xhistfeprem5,min=tplotran[0],max=tbinran[1])
yhistfeprem1=histogram(-strprepostfem1.tdiff/3600.,bin=3.,loc=xhistfeprem1,min=tplotran[0],max=tbinran[1])
yhistfeprec5=histogram(-strprepostfec5.tdiff/3600.,bin=3.,loc=xhistfeprec5,min=tplotran[0],max=tbinran[1])
yhistfeprec1=histogram(-strprepostfec1.tdiff/3600.,bin=3.,loc=xhistfeprec1,min=tplotran[0],max=tbinran[1])

;W.Time distributions
yhistfepostm5=histogram(strprepostfem5.ftdiff/3600.,bin=3.,loc=xhistfepostm5,min=tplotran[0],max=tbinran[1])
yhistfepostm1=histogram(strprepostfem1.ftdiff/3600.,bin=3.,loc=xhistfepostm1,min=tplotran[0],max=tbinran[1])
yhistfepostc5=histogram(strprepostfec5.ftdiff/3600.,bin=3.,loc=xhistfepostc5,min=tplotran[0],max=tbinran[1])
yhistfepostc1=histogram(strprepostfec1.ftdiff/3600.,bin=3.,loc=xhistfepostc1,min=tplotran[0],max=tbinran[1])

;Uncertainties of each distribution assuming poisson
uyhistfeprem5=sqrt(yhistfeprem5)
uyhistfepostm5=sqrt(yhistfepostm5)
uyhistfeprem1=sqrt(yhistfeprem1)
uyhistfepostm1=sqrt(yhistfepostm1)
uyhistfeprec5=sqrt(yhistfeprec5)
uyhistfepostc5=sqrt(yhistfepostc5)
uyhistfeprec1=sqrt(yhistfeprec1)
uyhistfepostc1=sqrt(yhistfepostc1)

;save structure with all arrays
xhistfepre=xhistfeprec1 & xhistfepost=xhistfepostc1
simstrfe={strprepostfec1:strprepostfec1,strprepostflc1:strprepostflc1, $
strprepostfec5:strprepostfec5,strprepostflc5:strprepostflc5, $
strprepostfem1:strprepostfem1,strprepostflm1:strprepostflm1, $
strprepostfem5:strprepostfem5,strprepostflm5:strprepostflm5, $
yhistfeprec1:yhistfeprec1,yhistfepostc1:yhistfepostc1, $
yhistfeprec5:yhistfeprec5,yhistfepostc5:yhistfepostc5, $
yhistfeprem1:yhistfeprem1,yhistfepostm1:yhistfepostm1, $
yhistfeprem5:yhistfeprem5,yhistfepostm5:yhistfepostm5, $
xhistfepre:xhistfepre, xhistfepost:xhistfepost}
save,simstrfe,file=datapath+'simulate_fe_fl_pair_relation-festruct.sav'

stop

;Plot >M5
!p.charsize=1.4
!p.multi=[0,1,2]
!p.background=255
!p.color=0

setplotenv,file=epspath+'simulate_fe-fl_pair_relation-compare_fl_fe_pre_post_gtm5_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostm5,yhistfepostm5,ps=10,yran=[0,20],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfeprem5,yhistfeprem5,ps=10,color=180
plot,xhistfepostm5,yhistfepostm5-yhistfeprem5,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-20,20],/ysty,xtit='M5-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistfepostm5, yhistfepostm5-yhistfeprem5, fltarr(n_elements(xhistfepostm5)), (uyhistfepostm5+uyhistfeprem5),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'simulate_fe-fl_pair_relation-compare_fl_fe_pre_post_gtm5_byfe'

stop

;Plot >M1
setplotenv,file=epspath+'simulate_fe-fl_pair_relation-compare_fl_fe_pre_post_gtm1_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostm1,yhistfepostm1,ps=10,yran=[0,70],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfeprem1,yhistfeprem1,ps=10,color=180
plot,xhistfepostm1,yhistfepostm1-yhistfeprem1,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-60,60],/ysty,xtit='M1-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistfepostm1, yhistfepostm1-yhistfeprem1, fltarr(n_elements(xhistfepostm1)), (uyhistfepostm1+uyhistfeprem1),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'simulate_fe-fl_pair_relation-compare_fl_fe_pre_post_gtm1_byfe'

stop

;Plot >C5
setplotenv,file=epspath+'simulate_fe-fl_pair_relation-compare_fl_fe_pre_post_gtc5_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostc5,yhistfepostc5,ps=10,yran=[0,120],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfeprec5,yhistfeprec5,ps=10,color=180
plot,xhistfepostc5,yhistfepostc5-yhistfeprec5,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-60,60],/ysty,xtit='C5-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistfepostc5, yhistfepostc5-yhistfeprec5, fltarr(n_elements(xhistfepostc5)), (uyhistfepostc5+uyhistfeprec5),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'simulate_fe-fl_pair_relation-compare_fl_fe_pre_post_gtc5_byfe'

stop

;Plot >C1
setplotenv,file=epspath+'simulate_fe-fl_pair_relation-compare_fl_fe_pre_post_gtc1_byfe.eps',/ps
!p.charsize=1.4
plot,xhistfepostc1,yhistfepostc1,ps=10,yran=[0,600],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfeprec1,yhistfeprec1,ps=10,color=180
plot,xhistfepostc1,yhistfepostc1-yhistfeprec1,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-150,150],/ysty,xtit='C1-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistfepostc1, yhistfepostc1-yhistfeprec1, fltarr(n_elements(xhistfepostc1)), (uyhistfepostc1+uyhistfeprec1),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'simulate_fe-fl_pair_relation-compare_fl_fe_pre_post_gtc1_byfe'

stop

;Plot the y-axis as the number of events in a given bin / total number of event pairs within 5 days
;Plot >M5, Fractional--------------------------------------------------------->
!p.charsize=2
!p.multi=[0,1,2]
!p.background=255
!p.color=0
!p.thick=5
!x.thick=2
!y.thick=2

stop

npostm5=float(n_elements(strprepostfem5.ftdiff))
npostm5hist=float(total(yhistfepostm5))
nprem5=float(n_elements(strprepostfem5.tdiff))
nprem5hist=float(total(yhistfeprem5))
yhistfeavgm5=(yhistfepostm5+yhistfeprem5)/2.

setplotenv,file=epspath+'simulate_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtm5_byfe.eps',/ps
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
;window_capture,file=plotpath+'simulate_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtm5_byfe'

stop

;Plot >M1, Fractional
npostm1=float(n_elements(strprepostfem1.ftdiff))
npostm1hist=float(total(yhistfepostm1))
nprem1=float(n_elements(strprepostfem1.tdiff))
nprem1hist=float(total(yhistfeprem1))
yhistfeavgm1=(yhistfepostm1+yhistfeprem1)/2.

setplotenv,file=epspath+'simulate_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtm1_byfe.eps',/ps,thick=5,xthick=7,ythick=7
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
;window_capture,file=plotpath+'simulate_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtm1_byfe'


stop

;Plot >C5, Fractional
npostc5=float(n_elements(strprepostfec5.ftdiff))
npostc5hist=float(total(yhistfepostc5))
nprec5=float(n_elements(strprepostfec5.tdiff))
nprec5hist=float(total(yhistfeprec5))
yhistfeavgc5=(yhistfepostc5+yhistfeprec5)/2.

setplotenv,file=epspath+'simulate_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtc5_byfe.eps',/ps
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
;window_capture,file=plotpath+'simulate_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtc5_byfe'

stop

;Plot >C1, Fractional
npostc1=float(n_elements(strprepostfec1.ftdiff))
npostc1hist=float(total(yhistfepostc1))
nprec1=float(n_elements(strprepostfec1.tdiff))
nprec1hist=float(total(yhistfeprec1))
yhistfeavgc1=(yhistfepostc1+yhistfeprec1)/2.

setplotenv,file=epspath+'simulate_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtc1_byfe.eps',/ps
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
;window_capture,file=plotpath+'simulate_fe-fl_pair_relation-frac_compare_fl_fe_pre_post_gtc1_byfe'


stop


stop






;Flare is the reference------------------------->

;W.Time distributions
yhistflprem5=histogram(-strprepostflm5.tdiff/3600.,bin=3.,loc=xhistflprem5,min=tplotran[0],max=tbinran[1])
yhistflprem1=histogram(-strprepostflm1.tdiff/3600.,bin=3.,loc=xhistflprem1,min=tplotran[0],max=tbinran[1])
yhistflprec5=histogram(-strprepostflc5.tdiff/3600.,bin=3.,loc=xhistflprec5,min=tplotran[0],max=tbinran[1])
yhistflprec1=histogram(-strprepostflc1.tdiff/3600.,bin=3.,loc=xhistflprec1,min=tplotran[0],max=tbinran[1])

;W.Time distributions
yhistflpostm5=histogram(strprepostflm5.ftdiff/3600.,bin=3.,loc=xhistflpostm5,min=tplotran[0],max=tbinran[1])
yhistflpostm1=histogram(strprepostflm1.ftdiff/3600.,bin=3.,loc=xhistflpostm1,min=tplotran[0],max=tbinran[1])
yhistflpostc5=histogram(strprepostflc5.ftdiff/3600.,bin=3.,loc=xhistflpostc5,min=tplotran[0],max=tbinran[1])
yhistflpostc1=histogram(strprepostflc1.ftdiff/3600.,bin=3.,loc=xhistflpostc1,min=tplotran[0],max=tbinran[1])

;Uncertainties of each distribution assuming poisson
uyhistflprem5=sqrt(yhistflprem5)
uyhistflpostm5=sqrt(yhistflpostm5)
uyhistflprem1=sqrt(yhistflprem1)
uyhistflpostm1=sqrt(yhistflpostm1)
uyhistflprec5=sqrt(yhistflprec5)
uyhistflpostc5=sqrt(yhistflpostc5)
uyhistflprec1=sqrt(yhistflprec1)
uyhistflpostc1=sqrt(yhistflpostc1)

;save structure with all arrays
xhistflpre=xhistflprec1 & xhistflpost=xhistflpostc1
simstrfl={strprepostflc1:strprepostflc1,strprepostfec1:strprepostfec1, $
strprepostflc5:strprepostflc5,strprepostfec5:strprepostfec5, $
strprepostflm1:strprepostflm1,strprepostfem1:strprepostfem1, $
strprepostflm5:strprepostflm5,strprepostfem5:strprepostfem5, $
yhistflprec1:yhistflprec1,yhistflpostc1:yhistflpostc1, $
yhistflprec5:yhistflprec5,yhistflpostc5:yhistflpostc5, $
yhistflprem1:yhistflprem1,yhistflpostm1:yhistflpostm1, $
yhistflprem5:yhistflprem5,yhistflpostm5:yhistflpostm5, $
xhistflpre:xhistflpre, xhistflpost:xhistflpost}
save,simstrfl,file=datapath+'simulate_fe_fl_pair_relation-flstruct.sav'

;Plot >M5

setplotenv,file=epspath+'simulate_fl-fe_pair_relation-compare_fl_fe_pre_post_gtm5_byfe.eps',/ps
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
;window_capture,file=plotpath+'simulate_fl-fe_pair_relation-compare_fl_fe_pre_post_gtm5_byfe'

stop

;Plot >M1

setplotenv,file=epspath+'simulate_fl-fe_pair_relation-compare_fl_fe_pre_post_gtm1_byfe.eps',/ps,thick=5,xthick=7,ythick=7
!p.charsize=1.4
!p.multi=[0,1,2]
plot,xhistflprem1,yhistflprem1,ps=10,yran=[0,25],xran=[0,24],ytit='# (FE. aft.=blk, bef.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistflpostm1,yhistflpostm1,ps=10,color=180
plot,xhistflpostm1,yhistflprem1-yhistflpostm1,xran=[0,24],ytit='(FE after Flare) - (FE before Flare)',ps=10,yran=[-20,20],/ysty,xtit='M1-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprem1, yhistflprem1-yhistflpostm1, fltarr(n_elements(xhistflprem1)), (uyhistflprem1+uyhistflpostm1),ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'simulate_fl-fe_pair_relation-compare_fl_fe_pre_post_gtm1_byfe'

stop

;Plot >C5

setplotenv,file=epspath+'simulate_fl-fe_pair_relation-compare_fl_fe_pre_post_gtc5_byfe.eps',/ps
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
;window_capture,file=plotpath+'simulate_fl-fe_pair_relation-compare_fl_fe_pre_post_gtc5_byfe'

stop

;Plot >C1

setplotenv,file=epspath+'simulate_fl-fe_pair_relation-compare_fl_fe_pre_post_gtc1_byfe.eps',/ps
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
;window_capture,file=plotpath+'simulate_fl-fe_pair_relation-compare_fl_fe_pre_post_gtc1_byfe'


stop

;Plot >M5; Fractional--------------------------------------------------------->

npostm5=float(n_elements(strprepostflm5.ftdiff))
npostm5hist=float(total(yhistflpostm5))
nprem5=float(n_elements(strprepostflm5.tdiff))
nprem5hist=float(total(yhistflprem5))
yhistflavgm5=(yhistflpostm5+yhistflprem5)/2.

setplotenv,file=epspath+'simulate_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtm5_byfe.eps',/ps
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
;window_capture,file=plotpath+'simulate_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtm5_byfe'

stop

;Plot >M1; Fractional

npostm1=float(n_elements(strprepostflm1.ftdiff))
npostm1hist=float(total(yhistflpostm1))
nprem1=float(n_elements(strprepostflm1.tdiff))
nprem1hist=float(total(yhistflprem1))
yhistflavgm1=(yhistflpostm1+yhistflprem1)/2.

setplotenv,file=epspath+'simulate_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtm1_byfe.eps',/ps,thick=5,xthick=7,ythick=7
!p.charsize=1.4
plot,xhistflprem1,yhistflprem1/nprem1,ps=10,yran=[0,.15],ytit='# / Ntot',ymarg=[2,1],xran=[0,24],xtit='',xtickname=strarr(10)+' ',ytickname=[' ']
sharpcorners
oploterr,xhistflprem1,yhistflprem1/nprem1,uyhistflprem1/nprem1,ps=4
legend,['FE Before','FE After'],lines=[0,0],color=[0,180],/top,/right
;vline,3,lines=2,color=150
oplot,xhistflpostm1,yhistflpostm1/npostm1,ps=10,color=180
!p.color=180
oploterr,xhistflpostm1,yhistflpostm1/npostm1,uyhistflpostm1/npostm1,ps=4
!p.color=0
;xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistflpostm1[1:7])/npostm1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflpostm1[8:*])/npostm1),color=180,/norm
;xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistflprem1[1:7])/nprem1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflprem1[8:*])/nprem1),color=0,/norm
plot,xhistflpostm1,(yhistflprem1-yhistflpostm1)/yhistflavgm1,ytit='(FE After - FE Before) / Avg(Before,After)',ps=10,yran=[-2,2],/ysty,xtit='M1-FE Wait-time [hrs]',xran=[0,24],ymarg=[5,-2]
sharpcorners
;vline,3,lines=2,color=150
hline,0,color=150
oploterr, xhistflprem1, (yhistflprem1-yhistflpostm1)/yhistflavgm1, fltarr(n_elements(xhistflprem1)), (uyhistflprem1+uyhistflpostm1)/yhistflavgm1,ps=4,ymarg=2
closeplotenv
;window_capture,file=plotpath+'simulate_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtm1_byfe'

stop

;Plot >C5; Fractional

npostc5=float(n_elements(strprepostflc5.ftdiff))
npostc5hist=float(total(yhistflpostc5))
nprec5=float(n_elements(strprepostflc5.tdiff))
nprec5hist=float(total(yhistflprec5))
yhistflavgc5=(yhistflpostc5+yhistflprec5)/2.

setplotenv,file=epspath+'simulate_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtc5_byfe.eps',/ps
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
;window_capture,file=plotpath+'simulate_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtc5_byfe'

stop

;Plot >C1; Fractional

npostc1=float(n_elements(strprepostflc1.ftdiff))
npostc1hist=float(total(yhistflpostc1))
nprec1=float(n_elements(strprepostflc1.tdiff))
nprec1hist=float(total(yhistflprec1))
yhistflavgc1=(yhistflpostc1+yhistflprec1)/2.

setplotenv,file=epspath+'simulate_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtc1_byfe.eps',/ps
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
;window_capture,file=plotpath+'simulate_fl-fe_pair_relation-frac_compare_fl_fe_pre_post_gtc1_byfe'

stop

;Make a plot summarising the comparison table
simsumarrfe=[[(yhistfepostc1-yhistfeprec1)/yhistfeavgc1],[(yhistfepostc5-yhistfeprec5)/yhistfeavgc5],[(yhistfepostm1-yhistfeprem1)/yhistfeavgm1],[(yhistfepostm5-yhistfeprem5)/yhistfeavgm5]]
simsumarrfl=[[(yhistflpostc1-yhistflprec1)/yhistflavgc1],[(yhistflpostc5-yhistflprec5)/yhistflavgc5],[(yhistflpostm1-yhistflprem1)/yhistflavgm1],[(yhistflpostm5-yhistflprem5)/yhistflavgm5]]
!p.multi=[0,1,2]
plot_image,simsumarrfe[0:9,*],scl=[0.333,1],ytit='FE ref'
plot_image,simsumarrfl[0:9,*],scl=[0.333,1],ytit='Flare ref'
save,simsumarrfe,simsumarrfl,file=datapath+'simulate_fe_fl_pair_relation-summaryarray'



;!!!! Re code the following to work for the above simulated time series
;!!! then make the corresponding plots, same as the observed ones, and see if they come out the same.

;!!! then run the boot strapping and do the same

;!!! then run the FE movies and pick out better positions... what do I observe? connections??


stop

end
