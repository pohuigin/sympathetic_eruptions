pro simulate_fe_fl_pair_bootstrap, res1tore=res1tore

;do boot strapping and reorder wait times, run comparisons FE-FL, and see if we get the same plots.

root='~/science/projects/stereo_sympathetic_flaring/'
datapath=root+'data/'
plotpath=root+'plots/determine_fe-fl_pair_relation/'
epspath=root+'plots/determine_fe-fl_pair_relation/eps/'
;epspath=root+'eps/'

restore,datapath+'determine_fe-fl_pair_relation.res1.sav',/ver

;Observed timeseries
tc1 = timle[wle[wc1]]
tc5 = timle[wle[wc5]]
tM1 = timle[wle[wM1]]
tM5 = timle[wle[wM5]]
tfe = timfe[wfe[wfe2]]

twc1 = flare_wait_time(tc1)
twc5 = flare_wait_time(tc5)
twM1 = flare_wait_time(tM1)
twm5 = flare_wait_time(tm5)
twfe = flare_wait_time(tfe)

nc1=n_elements(tc1)
nc5=n_elements(tc5)
nm1=n_elements(tm1)
nm5=n_elements(tm5)
nfe=n_elements(tfe)

nran=1000.
tplotran=[0.,24.] ;limit pairs to +-24
tbinran=[0.,7.*24.]

twranc1=(rebin(twc1,nc1,nran))
twranc5=(rebin(twc5,nc5,nran))
twranm1=(rebin(twm1,nm1,nran))
twranm5=(rebin(twm5,nm5,nran))
twranfe=(rebin(twfe,nfe,nran))

tranc1=fltarr(nc1,nran)
tranc5=fltarr(nc5,nran)
tranm1=fltarr(nm1,nran)
tranm5=fltarr(nm5,nran)
tranfe=fltarr(nfe,nran)

if not keyword_set(res1tore) then begin

for i=0,nran-1 do begin

	randexc1=sort(randomu(seed,nc1))
	randexc5=sort(randomu(seed,nc5))
	randexm1=sort(randomu(seed,nm1))
	randexm5=sort(randomu(seed,nm5))
	randexfe=sort(randomu(seed,nfe))

	twranc1[*,i]=(twranc1[*,i])[randexc1]
	twranc5[*,i]=(twranc5[*,i])[randexc5]
	twranm1[*,i]=(twranm1[*,i])[randexm1]
	twranm5[*,i]=(twranm5[*,i])[randexm5]
	twranfe[*,i]=(twranfe[*,i])[randexfe]

	tranc1[*,i]=flare_wait_time(twranc1[*,i],/inverse)
	tranc5[*,i]=flare_wait_time(twranc5[*,i],/inverse)
	tranm1[*,i]=flare_wait_time(twranm1[*,i],/inverse)
	tranm5[*,i]=flare_wait_time(twranm5[*,i],/inverse)
	tranfe[*,i]=flare_wait_time(twranfe[*,i],/inverse)

;Match up pairs of bootstrapped events
	istrprepostfec1=simulate_fe_fl_pair_relation_link_events_pre_post(tranfe[*,i], tranc1[*,i],/nopos)
	istrprepostfec5=simulate_fe_fl_pair_relation_link_events_pre_post(tranfe[*,i], tranc5[*,i],/nopos)
	istrprepostfem1=simulate_fe_fl_pair_relation_link_events_pre_post(tranfe[*,i], tranm1[*,i],/nopos)
	istrprepostfem5=simulate_fe_fl_pair_relation_link_events_pre_post(tranfe[*,i], tranm5[*,i],/nopos)

	istrprepostflc1=simulate_fe_fl_pair_relation_link_events_pre_post(tranc1[*,i],tranfe[*,i],/nopos)
	istrprepostflc5=simulate_fe_fl_pair_relation_link_events_pre_post(tranc5[*,i],tranfe[*,i], /nopos)
	istrprepostflm1=simulate_fe_fl_pair_relation_link_events_pre_post(tranm1[*,i],tranfe[*,i], /nopos)
	istrprepostflm5=simulate_fe_fl_pair_relation_link_events_pre_post(tranm5[*,i],tranfe[*,i], /nopos)

;W.Time distributions
;FE Reference
	iyhistfeprem5=histogram(-istrprepostfem5.tdiff/3600.,bin=3.,loc=xhistfeprem5,min=tplotran[0],max=tbinran[1])
	iyhistfeprem1=histogram(-istrprepostfem1.tdiff/3600.,bin=3.,loc=xhistfeprem1,min=tplotran[0],max=tbinran[1])
	iyhistfeprec5=histogram(-istrprepostfec5.tdiff/3600.,bin=3.,loc=xhistfeprec5,min=tplotran[0],max=tbinran[1])
	iyhistfeprec1=histogram(-istrprepostfec1.tdiff/3600.,bin=3.,loc=xhistfeprec1,min=tplotran[0],max=tbinran[1])

	iyhistfepostm5=histogram(istrprepostfem5.ftdiff/3600.,bin=3.,loc=xhistfepostm5,min=tplotran[0],max=tbinran[1])
	iyhistfepostm1=histogram(istrprepostfem1.ftdiff/3600.,bin=3.,loc=xhistfepostm1,min=tplotran[0],max=tbinran[1])
	iyhistfepostc5=histogram(istrprepostfec5.ftdiff/3600.,bin=3.,loc=xhistfepostc5,min=tplotran[0],max=tbinran[1])
	iyhistfepostc1=histogram(istrprepostfec1.ftdiff/3600.,bin=3.,loc=xhistfepostc1,min=tplotran[0],max=tbinran[1])

	if i eq 0 then begin
		yhistfeprem5arr=iyhistfeprem5
		yhistfeprem1arr=iyhistfeprem1
		yhistfeprec5arr=iyhistfeprec5
		yhistfeprec1arr=iyhistfeprec1
	
		yhistfepostm5arr=iyhistfepostm5
		yhistfepostm1arr=iyhistfepostm1
		yhistfepostc5arr=iyhistfepostc5
		yhistfepostc1arr=iyhistfepostc1
	endif else begin
		yhistfeprem5arr=[[yhistfeprem5arr],[iyhistfeprem5]]
		yhistfeprem1arr=[[yhistfeprem1arr],[iyhistfeprem1]]
		yhistfeprec5arr=[[yhistfeprec5arr],[iyhistfeprec5]]
		yhistfeprec1arr=[[yhistfeprec1arr],[iyhistfeprec1]]

		yhistfepostm5arr=[[yhistfepostm5arr],[iyhistfepostm5]]
		yhistfepostm1arr=[[yhistfepostm1arr],[iyhistfepostm1]]
		yhistfepostc5arr=[[yhistfepostc5arr],[iyhistfepostc5]]
		yhistfepostc1arr=[[yhistfepostc1arr],[iyhistfepostc1]]
	endelse

;FL Reference
	iyhistflprem5=histogram(-istrprepostflm5.tdiff/3600.,bin=3.,loc=xhistflprem5,min=tplotran[0],max=tbinran[1])
	iyhistflprem1=histogram(-istrprepostflm1.tdiff/3600.,bin=3.,loc=xhistflprem1,min=tplotran[0],max=tbinran[1])
	iyhistflprec5=histogram(-istrprepostflc5.tdiff/3600.,bin=3.,loc=xhistflprec5,min=tplotran[0],max=tbinran[1])
	iyhistflprec1=histogram(-istrprepostflc1.tdiff/3600.,bin=3.,loc=xhistflprec1,min=tplotran[0],max=tbinran[1])

	iyhistflpostm5=histogram(istrprepostflm5.ftdiff/3600.,bin=3.,loc=xhistflpostm5,min=tplotran[0],max=tbinran[1])
	iyhistflpostm1=histogram(istrprepostflm1.ftdiff/3600.,bin=3.,loc=xhistflpostm1,min=tplotran[0],max=tbinran[1])
	iyhistflpostc5=histogram(istrprepostflc5.ftdiff/3600.,bin=3.,loc=xhistflpostc5,min=tplotran[0],max=tbinran[1])
	iyhistflpostc1=histogram(istrprepostflc1.ftdiff/3600.,bin=3.,loc=xhistflpostc1,min=tplotran[0],max=tbinran[1])

	if i eq 0 then begin
		yhistflprem5arr=iyhistflprem5
		yhistflprem1arr=iyhistflprem1
		yhistflprec5arr=iyhistflprec5
		yhistflprec1arr=iyhistflprec1
	
		yhistflpostm5arr=iyhistflpostm5
		yhistflpostm1arr=iyhistflpostm1
		yhistflpostc5arr=iyhistflpostc5
		yhistflpostc1arr=iyhistflpostc1
	endif else begin
		yhistflprem5arr=[[yhistflprem5arr],[iyhistflprem5]]
		yhistflprem1arr=[[yhistflprem1arr],[iyhistflprem1]]
		yhistflprec5arr=[[yhistflprec5arr],[iyhistflprec5]]
		yhistflprec1arr=[[yhistflprec1arr],[iyhistflprec1]]

		yhistflpostm5arr=[[yhistflpostm5arr],[iyhistflpostm5]]
		yhistflpostm1arr=[[yhistflpostm1arr],[iyhistflpostm1]]
		yhistflpostc5arr=[[yhistflpostc5arr],[iyhistflpostc5]]
		yhistflpostc1arr=[[yhistflpostc1arr],[iyhistflpostc1]]
	endelse

endfor

xhistfl=xhistflprem5
xhistfe=xhistfeprem5

save,yhistfeprem5arr,yhistfeprem1arr,yhistfeprec5arr,yhistfeprec1arr, $
	yhistfepostm5arr,yhistfepostm1arr,yhistfepostc5arr,yhistfepostc1arr, $
	yhistflprem5arr,yhistflprem1arr,yhistflprec5arr,yhistflprec1arr, $
	yhistflpostm5arr,yhistflpostm1arr,yhistflpostc5arr,yhistflpostc1arr, $
	xhistfl,xhistfe, $
	file=datapath+'simulate_fe_fl_pair_bootstrap-res1.sav'
	
endif else restore,datapath+'simulate_fe_fl_pair_bootstrap-res1.sav',/ver

bootstr={yhistfeprem5arr:yhistfeprem5arr,yhistfeprem1arr:yhistfeprem1arr,yhistfeprec5arr:yhistfeprec5arr,yhistfeprec1arr:yhistfeprec1arr, $
	yhistfepostm5arr:yhistfepostm5arr,yhistfepostm1arr:yhistfepostm1arr,yhistfepostc5arr:yhistfepostc5arr,yhistfepostc1arr:yhistfepostc1arr, $
	yhistflprem5arr:yhistflprem5arr,yhistflprem1arr:yhistflprem1arr,yhistflprec5arr:yhistflprec5arr,yhistflprec1arr:yhistflprec1arr, $
	yhistflpostm5arr:yhistflpostm5arr,yhistflpostm1arr:yhistflpostm1arr,yhistflpostc5arr:yhistflpostc5arr,yhistflpostc1arr:yhistflpostc1arr, $
	xhistfl:xhistfl,xhistfe:xhistfe}
save,bootstr,file=datapath+'simulate_fe_fl_pair_bootstrap-struct.sav'


!p.multi=[0,2,2]
loadct,0,/sil
setcolors,/sys,/sil
!p.color=0
!p.background=255

plot_hist,(tc1-min(tc1))/3600./24.,charsize=1.4
plot_hist,tranc1[*,0]/3600./24.,/oplot,color=!red
plot_hist,(tc5-min(tc5))/3600./24.,charsize=1.4
plot_hist,tranc5[*,0]/3600./24.,/oplot,color=!red

plot_hist,(tm1-min(tm1))/3600./24.,charsize=1.4
plot_hist,tranm1[*,0]/3600./24.,/oplot,color=!red
plot_hist,(tm5-min(tm5))/3600./24.,charsize=1.4
plot_hist,tranm1[*,0]/3600./24.,/oplot,color=!red

;plot_hist,(tfe-min(tfe))/3600./24.,/oplot,color=!gray
;plot_hist,tranfe[*,0]/3600./24.,/oplot,color=!blue

stop


;Plot the difference between each pair---------------------------------------->

tplotran=[0.,24.] ;limit pairs to +-24
tbinran=[0.,7.*24.]

;Filament Eruption is the reference------------->

!p.multi=[0,1,2]
!p.background=255
!p.color=0

;Plot >M5
yhistfepostm5=total(yhistfepostm5arr,2)/nran
yhistfeprem5=total(yhistfeprem5arr,2)/nran

setplotenv,file=epspath+'simulate_fe-fl_pair_bootstrap-compare_fl_fe_pre_post_gtm5_byfe.eps',/ps,ys=7,xs=7
!p.charsize=1.4
plot,xhistfe,yhistfepostm5,ps=10,yran=[0,10],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfe,yhistfeprem5,ps=10,color=180
plot,xhistfe,yhistfepostm5-yhistfeprem5,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-1,1],/ysty,xtit='M5-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
;oploterr, xhistfe, yhistfepostm5-yhistfeprem5, fltarr(n_elements(xhistfepostm5)), sqrt(uyhistfepostm5^2.+uyhistfeprem5^2.),ps=4,ymarg=2
;window_capture,file=plotpath+'simulate_fe-fl_pair_bootstrap-compare_fl_fe_pre_post_gtm5_byfe'
closeplotenv

stop

;Plot >M1
yhistfepostm1=total(yhistfepostm1arr,2)/nran
yhistfeprem1=total(yhistfeprem1arr,2)/nran
;uncertainties
uyhistfeprem1arr=fltarr(n_elements(yhistfeprem1arr[*,0]))
uyhistfepostm1arr=fltarr(n_elements(yhistfepostm1arr[*,0]))
for i=0,n_elements(yhistfeprem1arr[*,0])-1 do begin
	uyhistfeprem1arr[i]=stddev(yhistfeprem1arr[i,*])
	uyhistfepostm1arr[i]=stddev(yhistfepostm1arr[i,*])
endfor


setplotenv,file=epspath+'simulate_fe-fl_pair_bootstrap-compare_fl_fe_pre_post_gtm1_byfe.eps',/ps,thick=5,xthick=7,ythick=7
!p.charsize=1.4
plot,xhistfe,yhistfepostm1,ps=10,yran=[0,20],ytit='# Pairs',ymarg=[2,1],xran=[0,24],xtit='',xtickname=strarr(10)+' ',ytickname=[' ']
oploterr,xhistfe,yhistfepostm1,uyhistfepostm1arr,ps=4
legend,['Flare Before','Flare After'],lines=[0,0],color=[0,180],/bottom,/left
;vline,3,lines=2,color=150
oplot,xhistfe,yhistfeprem1,ps=10,color=180
oploterr,xhistfe,yhistfeprem1,uyhistfeprem1arr,ps=4,color=180
plot,xhistfe,yhistfepostm1-yhistfeprem1,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-10,10],/ysty,xtit='M1-FE Wait-time [hrs]',ymarg=[5,-2]
oploterr,xhistfe,yhistfepostm1-yhistfeprem1,uyhistfeprem1arr+uyhistfepostm1arr,ps=4
;vline,3,lines=2,color=150
hline,0,thick=2,color=150
;oploterr, xhistfepostm1, yhistfepostm1-yhistfeprem1, fltarr(n_elements(xhistfepostm1)), sqrt(uyhistfepostm1^2.+uyhistfeprem1^2.),ps=4,ymarg=2
;window_capture,file=plotpath+'simulate_fe-fl_pair_bootstrap-compare_fl_fe_pre_post_gtm1_byfe'
closeplotenv

stop

;Plot >C5
yhistfepostc5=total(yhistfepostc5arr,2)/nran
yhistfeprec5=total(yhistfeprec5arr,2)/nran

setplotenv,file=epspath+'simulate_fe-fl_pair_bootstrap-compare_fl_fe_pre_post_gtc5_byfe.eps',/ps,ys=7,xs=7
!p.charsize=1.4
plot,xhistfe,yhistfepostc5,ps=10,yran=[0,20],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfe,yhistfeprec5,ps=10,color=180
plot,xhistfe,yhistfepostc5-yhistfeprec5,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-1,1],/ysty,xtit='C5-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
;oploterr, xhistfepostc5, yhistfepostc5-yhistfeprec5, fltarr(n_elements(xhistfepostc5)), sqrt(uyhistfepostc5^2.+uyhistfeprec5^2.),ps=4,ymarg=2
;window_capture,file=plotpath+'simulate_fe-fl_pair_bootstrap-compare_fl_fe_pre_post_gtc5_byfe'
closeplotenv

stop

;Plot >C1
yhistfepostc1=total(yhistfepostc1arr,2)/nran
yhistfeprec1=total(yhistfeprec1arr,2)/nran

setplotenv,file=epspath+'simulate_fe-fl_pair_bootstrap-compare_fl_fe_pre_post_gtc1_byfe.eps',/ps,ys=7,xs=7
!p.charsize=1.4
plot,xhistfe,yhistfepostc1,ps=10,yran=[0,20],xran=[0,24],ytit='# (Fl. bef.=blk, aft.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfe,yhistfeprec1,ps=10,color=180
plot,xhistfe,yhistfepostc1-yhistfeprec5,xran=[0,24],ytit='(Flare before FE) - (Flare after FE)',ps=10,yran=[-1,1],/ysty,xtit='C1-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
;oploterr, xhistfepostc1, yhistfepostc1-yhistfeprec1, fltarr(n_elements(xhistfepostc1)), sqrt(uyhistfepostc1^2.+uyhistfeprec1^2.),ps=4,ymarg=2
;window_capture,file=plotpath+'simulate_fe-fl_pair_bootstrap-compare_fl_fe_pre_post_gtc1_byfe'
closeplotenv

stop

;-------------------------->
;With the boot strapping, we are averaging, so y-axis isn't the actual number of events anyhow
goto,skipfracplots


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

plot,xhistfepostm5,yhistfepostm5/npostm5,ps=10,yran=[0,.02],ytit='#/Ntot (Fl. bef.=blk, aft.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=180
oplot,xhistfeprem5,yhistfeprem5/nprem5,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistfepostm5[1:7])/npostm5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfepostm5[8:*])/npostm5),color=0,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistfeprem5[1:7])/nprem5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfeprem5[8:*])/nprem5),color=180,/norm
plot,xhistfepostm5,(yhistfepostm5-yhistfeprem5)/yhistfeavgm5,xran=[0,24],ytit='(Flare bef. - Flare aft.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='M5-FE Wait-time [hrs]'
vline,3,lines=2,color=180
hline,0,thick=2,color=150
oploterr, xhistfepostm5, (yhistfepostm5-yhistfeprem5)/yhistfeavgm5, fltarr(n_elements(xhistfepostm5)), sqrt(uyhistfepostm5^2.+uyhistfeprem5^2.)/yhistfeavgm5,ps=4,ymarg=2
window_capture,file=plotpath+'simulate_fe-fl_pair_bootstrap-frac_compare_fl_fe_pre_post_gtm5_byfe'

stop

;Plot >M1, Fractional
npostm1=float(n_elements(strprepostfem1.ftdiff))
npostm1hist=float(total(yhistfepostm1))
nprem1=float(n_elements(strprepostfem1.tdiff))
nprem1hist=float(total(yhistfeprem1))
yhistfeavgm1=(yhistfepostm1+yhistfeprem1)/2.

plot,xhistfepostm1,yhistfepostm1/npostm1,ps=10,yran=[0,.04],ytit='#/Ntot (Fl. bef.=blk, aft.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=180
oplot,xhistfeprem1,yhistfeprem1/nprem1,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistfepostm1[1:7])/npostm1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfepostm1[8:*])/npostm1),color=0,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistfeprem1[1:7])/nprem1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfeprem1[8:*])/nprem1),color=180,/norm
plot,xhistfepostm1,(yhistfepostm1-yhistfeprem1)/yhistfeavgm1,ytit='(Flare bef. - Flare aft.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='M1-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=180
hline,0,thick=2,color=150
oploterr, xhistfepostm1, (yhistfepostm1-yhistfeprem1)/yhistfeavgm1, fltarr(n_elements(xhistfepostm1)), sqrt(uyhistfepostm1^2.+uyhistfeprem1^2.)/yhistfeavgm1,ps=4,ymarg=2
window_capture,file=plotpath+'simulate_fe-fl_pair_bootstrap-frac_compare_fl_fe_pre_post_gtm1_byfe'


stop

;Plot >C5, Fractional
npostc5=float(n_elements(strprepostfec5.ftdiff))
npostc5hist=float(total(yhistfepostc5))
nprec5=float(n_elements(strprepostfec5.tdiff))
nprec5hist=float(total(yhistfeprec5))
yhistfeavgc5=(yhistfepostc5+yhistfeprec5)/2.

plot,xhistfepostc5,yhistfepostc5/npostc5,ps=10,yran=[0,.08],ytit='#/Ntot (Fl. bef.=blk, aft.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=180
oplot,xhistfeprec5,yhistfeprec5/nprec5,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistfepostc5[1:7])/npostc5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfepostc5[8:*])/npostc5),color=0,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistfeprec5[1:7])/nprec5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfeprec5[8:*])/nprec5),color=180,/norm
plot,xhistfepostc5,(yhistfepostc5-yhistfeprec5)/yhistfeavgc5,ytit='(Flare bef. - Flare aft.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='M5-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=180
hline,0,thick=2,color=150
oploterr, xhistfepostc5, (yhistfepostc5-yhistfeprec5)/yhistfeavgc5, fltarr(n_elements(xhistfepostc5)), sqrt(uyhistfepostc5^2.+uyhistfeprec5^2.)/yhistfeavgc5,ps=4,ymarg=2
window_capture,file=plotpath+'simulate_fe-fl_pair_bootstrap-frac_compare_fl_fe_pre_post_gtc5_byfe'

stop

;Plot >C1, Fractional
npostc1=float(n_elements(strprepostfec1.ftdiff))
npostc1hist=float(total(yhistfepostc1))
nprec1=float(n_elements(strprepostfec1.tdiff))
nprec1hist=float(total(yhistfeprec1))
yhistfeavgc1=(yhistfepostc1+yhistfeprec1)/2.

plot,xhistfepostc1,yhistfepostc1/npostc1,ps=10,yran=[0,.16],ytit='#/Ntot (Fl. bef.=blk, aft.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=180
oplot,xhistfeprec1,yhistfeprec1/nprec1,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistfepostc1[1:7])/npostc1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfepostc1[8:*])/npostc1),color=0,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistfeprec1[1:7])/nprec1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistfeprec1[8:*])/nprec1),color=180,/norm
plot,xhistfepostc1,(yhistfepostc1-yhistfeprec1)/yhistfeavgc1,ytit='(Flare bef. - Flare aft.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='M5-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=180
hline,0,thick=2,color=150
oploterr, xhistfepostc1, (yhistfepostc1-yhistfeprec1)/yhistfeavgc1, fltarr(n_elements(xhistfepostc1)), sqrt(uyhistfepostc1^2.+uyhistfeprec1^2.)/yhistfeavgc1,ps=4,ymarg=2
window_capture,file=plotpath+'simulate_fe-fl_pair_bootstrap-frac_compare_fl_fe_pre_post_gtc1_byfe'


stop

skipfracplots:













;Flare is the reference------------------------->

;Plot >M5

yhistflpostm5=total(yhistflpostm5arr,2)/nran
yhistflprem5=total(yhistflprem5arr,2)/nran

setplotenv,file=epspath+'simulate_fl-fe_pair_bootstrap-compare_fl_fe_pre_post_gtm5_byfe.eps',/ps,ys=7,xs=7
!p.charsize=1.4
plot,xhistfl,yhistflprem5,ps=10,yran=[0,10],xran=[0,24],ytit='# (FE. aft.=blk, bef.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfl,yhistflpostm5,ps=10,color=180
plot,xhistfl,yhistflprem5-yhistflpostm5,xran=[0,24],ytit='(FE after Flare) - (FE before Flare)',ps=10,yran=[-1,1],/ysty,xtit='M5-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
;oploterr, xhistflprem5, yhistflprem5-yhistflpostm5, fltarr(n_elements(xhistflprem5)), sqrt(uyhistflprem5^2.+uyhistflpostm5^2.),ps=4,ymarg=2
;window_capture,file=plotpath+'simulate_fl-fe_pair_bootstrap-compare_fl_fe_pre_post_gtm5_byfe'
closeplotenv

stop

;Plot >M1

yhistflpostm1=total(yhistflpostm1arr,2)/nran
yhistflprem1=total(yhistflprem1arr,2)/nran
;uncertainties
uyhistflprem1arr=fltarr(n_elements(yhistflprem1arr[*,0]))
uyhistflpostm1arr=fltarr(n_elements(yhistflpostm1arr[*,0]))
for i=0,n_elements(yhistflprem1arr[*,0])-1 do begin
	uyhistflprem1arr[i]=stddev(yhistflprem1arr[i,*])
	uyhistflpostm1arr[i]=stddev(yhistflpostm1arr[i,*])
endfor

setplotenv,file=epspath+'simulate_fl-fe_pair_bootstrap-compare_fl_fe_pre_post_gtm1_byfe.eps',/ps,thick=5,xthick=7,ythick=7
!p.charsize=1.4
plot,xhistfl,yhistflprem1,ps=10,yran=[0,20],xran=[0,24],ytit='# Pairs',ymarg=[2,1],xtit='',xtickname=strarr(10)+' ',ytickname=[' ']
oploterr,xhistfl,yhistflprem1,uyhistfeprem1arr,ps=4
legend,['Flare Before','Flare After'],lines=[0,0],color=[0,180],/bottom,/left
;vline,3,lines=2,color=150
oplot,xhistfl,yhistflpostm1,ps=10,color=180
plot,xhistfl,yhistflprem1-yhistflpostm1,xran=[0,24],ytit='(FE after Flare) - (FE before Flare)',ps=10,yran=[-10,10],/ysty,xtit='M1-FE Wait-time [hrs]',ymarg=[5,-2]
oploterr,xhistfl,yhistflprem1-yhistflpostm1,uyhistflprem1arr+uyhistflpostm1arr,ps=4
;vline,3,lines=2,color=150
hline,0,thick=2,color=150
;oploterr, xhistflprem1, yhistflprem1-yhistflpostm1, fltarr(n_elements(xhistflprem1)), sqrt(uyhistflprem1^2.+uyhistflpostm1^2.),ps=4,ymarg=2
;window_capture,file=plotpath+'simulate_fl-fe_pair_bootstrap-compare_fl_fe_pre_post_gtm1_byfe'
closeplotenv

stop

;Plot >C5

yhistflpostc5=total(yhistflpostc5arr,2)/nran
yhistflprec5=total(yhistflprec5arr,2)/nran

setplotenv,file=epspath+'simulate_fl-fe_pair_bootstrap-compare_fl_fe_pre_post_gtc5_byfe.eps',/ps,ys=7,xs=7
!p.charsize=1.4
plot,xhistfl,yhistflprec5,ps=10,yran=[0,10],xran=[0,24],ytit='# (FE. aft.=blk, bef.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfl,yhistflpostc5,ps=10,color=180
plot,xhistfl,yhistflprec5-yhistflpostc5,xran=[0,24],ytit='(FE after Flare) - (FE before Flare)',ps=10,yran=[-1,1],/ysty,xtit='C5-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
;oploterr, xhistflprec5, yhistflprec5-yhistflpostc5, fltarr(n_elements(xhistflprec5)), sqrt(uyhistflprec5^2.+uyhistflpostc5^2.),ps=4,ymarg=2
;window_capture,file=plotpath+'simulate_fl-fe_pair_bootstrap-compare_fl_fe_pre_post_gtc5_byfe'
closeplotenv

stop

;Plot >C1

yhistflpostc1=total(yhistflpostc1arr,2)/nran
yhistflprec1=total(yhistflprec1arr,2)/nran

setplotenv,file=epspath+'simulate_fl-fe_pair_bootstrap-compare_fl_fe_pre_post_gtc1_byfe.eps',/ps,ys=7,xs=7
!p.charsize=1.4
plot,xhistfl,yhistflprec1,ps=10,yran=[0,10],xran=[0,24],ytit='# (FE. aft.=blk, bef.=gray)',ymarg=1
vline,3,lines=2,color=150
oplot,xhistfl,yhistflpostc1,ps=10,color=180
plot,xhistfl,yhistflprec1-yhistflpostc1,xran=[0,24],ytit='(FE after Flare) - (FE before Flare)',ps=10,yran=[-1,1],/ysty,xtit='C1-FE Wait-time [hrs]'
vline,3,lines=2,color=150
hline,0,thick=2,color=150
;oploterr, xhistflprec1, yhistflprec1-yhistflpostc1, fltarr(n_elements(xhistflprec1)), sqrt(uyhistflprec1^2.+uyhistflpostc1^2.),ps=4,ymarg=2
;window_capture,file=plotpath+'simulate_fl-fe_pair_bootstrap-compare_fl_fe_pre_post_gtc1_byfe'
closeplotenv

stop



goto,skipfractionalplots2

;Plot >M5; Fractional--------------------------------------------------------->

npostm5=float(n_elements(strprepostflm5.ftdiff))
npostm5hist=float(total(yhistflpostm5))
nprem5=float(n_elements(strprepostflm5.tdiff))
nprem5hist=float(total(yhistflprem5))
yhistflavgm5=(yhistflpostm5+yhistflprem5)/2.

plot,xhistflprem5,yhistflprem5/nprem5,ps=10,yran=[0,.2],ytit='#/Ntot (FE. aft.=blk, bef.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=150
oplot,xhistflpostm5,yhistflpostm5/npostm5,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistflpostm5[1:7])/npostm5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflpostm5[8:*])/npostm5),color=180,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistflprem5[1:7])/nprem5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflprem5[8:*])/nprem5),color=0,/norm
plot,xhistflpostm5,(yhistflprem5-yhistflpostm5)/yhistflavgm5,ytit='(FE aft. - FE bef.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='M5-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprem5, (yhistflprem5-yhistflpostm5)/yhistflavgm5, fltarr(n_elements(xhistflprem5)), sqrt(uyhistflprem5^2.+uyhistflpostm5^2.)/yhistflavgm5,ps=4,ymarg=2
window_capture,file=plotpath+'simulate_fl-fe_pair_bootstrap-frac_compare_fl_fe_pre_post_gtm5_byfe'

stop

;Plot >M1; Fractional

npostm1=float(n_elements(strprepostflm1.ftdiff))
npostm1hist=float(total(yhistflpostm1))
nprem1=float(n_elements(strprepostflm1.tdiff))
nprem1hist=float(total(yhistflprem1))
yhistflavgm1=(yhistflpostm1+yhistflprem1)/2.

plot,xhistflprem1,yhistflprem1/nprem1,ps=10,yran=[0,.1],ytit='#/Ntot (FE. aft.=blk, bef.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=150
oplot,xhistflpostm1,yhistflpostm1/npostm1,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistflpostm1[1:7])/npostm1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflpostm1[8:*])/npostm1),color=180,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistflprem1[1:7])/nprem1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflprem1[8:*])/nprem1),color=0,/norm
plot,xhistflpostm1,(yhistflprem1-yhistflpostm1)/yhistflavgm1,ytit='(FE aft. - FE bef.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='M1-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprem1, (yhistflprem1-yhistflpostm1)/yhistflavgm1, fltarr(n_elements(xhistflprem1)), sqrt(uyhistflprem1^2.+uyhistflpostm1^2.)/yhistflavgm1,ps=4,ymarg=2
window_capture,file=plotpath+'simulate_fl-fe_pair_bootstrap-frac_compare_fl_fe_pre_post_gtm1_byfe'

stop

;Plot >C5; Fractional

npostc5=float(n_elements(strprepostflc5.ftdiff))
npostc5hist=float(total(yhistflpostc5))
nprec5=float(n_elements(strprepostflc5.tdiff))
nprec5hist=float(total(yhistflprec5))
yhistflavgc5=(yhistflpostc5+yhistflprec5)/2.

plot,xhistflprec5,yhistflprec5/nprec5,ps=10,yran=[0,.1],ytit='#/Ntot (FE. aft.=blk, bef.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=150
oplot,xhistflpostc5,yhistflpostc5/npostc5,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistflpostc5[1:7])/npostc5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflpostc5[8:*])/npostc5),color=180,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistflprec5[1:7])/nprec5)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflprec5[8:*])/nprec5),color=0,/norm
plot,xhistflpostc5,(yhistflprec5-yhistflpostc5)/yhistflavgc5,ytit='(FE aft. - FE bef.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='C5-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprec5, (yhistflprec5-yhistflpostc5)/yhistflavgc5, fltarr(n_elements(xhistflprec5)), sqrt(uyhistflprec5^2.+uyhistflpostc5^2.)/yhistflavgc5,ps=4,ymarg=2
window_capture,file=plotpath+'simulate_fl-fe_pair_bootstrap-frac_compare_fl_fe_pre_post_gtc5_byfe'

stop

;Plot >C1; Fractional

npostc1=float(n_elements(strprepostflc1.ftdiff))
npostc1hist=float(total(yhistflpostc1))
nprec1=float(n_elements(strprepostflc1.tdiff))
nprec1hist=float(total(yhistflprec1))
yhistflavgc1=(yhistflpostc1+yhistflprec1)/2.

plot,xhistflprec1,yhistflprec1/nprec1,ps=10,yran=[0,.1],ytit='#/Ntot (FE. aft.=blk, bef.=gray)',ymarg=1,xran=[0,24]
vline,3,lines=2,color=150
oplot,xhistflpostc1,yhistflpostc1/npostc1,ps=10,color=180
xyouts,0.15,0.9,'# 3->24hr = '+strtrim(total(yhistflpostc1[1:7])/npostc1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflpostc1[8:*])/npostc1),color=180,/norm
xyouts,0.15,0.85,'# 3->24hr = '+strtrim(total(yhistflprec1[1:7])/nprec1)+'; #/Ntot 24->168hr = '+strtrim(total(yhistflprec1[8:*])/nprec1),color=0,/norm
plot,xhistflpostc1,(yhistflprec1-yhistflpostc1)/yhistflavgc1,ytit='(FE aft. - FE bef.)/Avg(Bef.,Aft.)',ps=10,yran=[-2,2],/ysty,xtit='C5-FE Wait-time [hrs]',xran=[0,24]
vline,3,lines=2,color=150
hline,0,thick=2,color=150
oploterr, xhistflprec1, (yhistflprec1-yhistflpostc1)/yhistflavgc1, fltarr(n_elements(xhistflprec1)), sqrt(uyhistflprec1^2.+uyhistflpostc1^2.)/yhistflavgc1,ps=4,ymarg=2
window_capture,file=plotpath+'simulate_fl-fe_pair_bootstrap-frac_compare_fl_fe_pre_post_gtc1_byfe'


skipfractionalplots2:

stop




;Put plots from poisson stat sim and boot strap sim into pres. -> send to karel


;!!! then run the FE movies and pick out better positions... what do I observe? connections??







end
