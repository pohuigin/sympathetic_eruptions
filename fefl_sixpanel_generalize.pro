;Calculate KS for bootstrap for 4 flare thresholds.

pro fefl_sixpanel_generalize_bootks,bootstr,ksbootstr,ksbootstrarr,xran=xran

nhist=n_elements((bootstr.yhistfeprem5arr)[0,*])

ksbootstrarr=replicate(ksbootstr,nhist)

for i=0,nhist-1 do begin

;Bootstrapped 0-24 hours
kstwo, ((bootstr.yhistfeprem5arr)[*,i])[xran[0]:xran[1]], ((bootstr.yhistfepostm5arr)[*,i])[xran[0]:xran[1]], dfem5, pfem5
ksbootstrarr[i].dfem5=dfem5 & ksbootstrarr[i].pfem5=pfem5
kstwo, ((bootstr.yhistfeprem1arr)[*,i])[xran[0]:xran[1]], ((bootstr.yhistfepostm1arr)[*,i])[xran[0]:xran[1]], dfem1, pfem1
ksbootstrarr[i].dfem1=dfem1 & ksbootstrarr[i].pfem1=pfem1
kstwo, ((bootstr.yhistfeprec5arr)[*,i])[xran[0]:xran[1]], ((bootstr.yhistfepostc5arr)[*,i])[xran[0]:xran[1]], dfec5, pfec5
ksbootstrarr[i].dfec5=dfec5 & ksbootstrarr[i].pfec5=pfec5
kstwo, ((bootstr.yhistfeprec1arr)[*,i])[xran[0]:xran[1]], ((bootstr.yhistfepostc1arr)[*,i])[xran[0]:xran[1]], dfec1, pfec1
ksbootstrarr[i].dfec1=dfec1 & ksbootstrarr[i].pfec1=pfec1

kstwo, ((bootstr.yhistflprem5arr)[*,i])[xran[0]:xran[1]], ((bootstr.yhistflpostm5arr)[*,i])[xran[0]:xran[1]], dflm5, pflm5
ksbootstrarr[i].dflm5=dflm5 & ksbootstrarr[i].pflm5=pflm5
kstwo, ((bootstr.yhistflprem1arr)[*,i])[xran[0]:xran[1]], ((bootstr.yhistflpostm1arr)[*,i])[xran[0]:xran[1]], dflm1, pflm1
ksbootstrarr[i].dflm1=dflm1 & ksbootstrarr[i].pflm1=pflm1
kstwo, ((bootstr.yhistflprec5arr)[*,i])[xran[0]:xran[1]], ((bootstr.yhistflpostc5arr)[*,i])[xran[0]:xran[1]], dflc5, pflc5
ksbootstrarr[i].dflc5=dflc5 & ksbootstrarr[i].pflc5=pflc5
kstwo, ((bootstr.yhistflprec1arr)[*,i])[xran[0]:xran[1]], ((bootstr.yhistflpostc1arr)[*,i])[xran[0]:xran[1]], dflc1, pflc1
ksbootstrarr[i].dflc1=dflc1 & ksbootstrarr[i].pflc1=pflc1

endfor

ksbootstr.dfem5=mean(ksbootstrarr.dfem5) & ksbootstr.pfem5=mean(ksbootstrarr.pfem5)

ksbootstr.dfem1=mean(ksbootstrarr.dfem1) & ksbootstr.pfem1=mean(ksbootstrarr.pfem1)

ksbootstr.dfec5=mean(ksbootstrarr.dfec5) & ksbootstr.pfec5=mean(ksbootstrarr.pfec5)

ksbootstr.dfec1=mean(ksbootstrarr.dfec1) & ksbootstr.pfec1=mean(ksbootstrarr.pfec1)

ksbootstr.dflm5=mean(ksbootstrarr.dflm5) & ksbootstr.pflm5=mean(ksbootstrarr.pflm5)

ksbootstr.dflm1=mean(ksbootstrarr.dflm1) & ksbootstr.pflm1=mean(ksbootstrarr.pflm1)

ksbootstr.dflc5=mean(ksbootstrarr.dflc5) & ksbootstr.pflc5=mean(ksbootstrarr.pflc5)

ksbootstr.dflc1=mean(ksbootstrarr.dflc1) & ksbootstr.pflc1=mean(ksbootstrarr.pflc1)



end

;----------------------------------------------------------------------------->

;Make 6-panel plot with obs/sim/boot for 4 flare thresholds

pro fefl_sixpanel_plotdist,param,obsstrfe,obsstrfl,simstrfe,simstrfl,bootstr

cls=param.cls

!p.charsize=2
!p.multi=[0,2,6]
!p.background=255
!p.color=0
!p.charsize=1.4

setplotenv,file=param.epspath+'fefl_sixpanel'+cls+'.eps',/ps,xs=10,ys=15

;Black=Matched Event After
;Grey=Matched Event Before

;PRE/POST corresponds to the Matched event coming BEFORE/AFTER the Reference event

;N_BEFORE,AFTER
;Tau_REFERENCE,MATCHED

;------------------------------------------------------------------->

status=execute('npostcls=float(n_elements(obsstrfe.strprepostfe'+cls+'.ftdiff))')
status=execute('npostclshist=float(total(obsstrfe.yhistfepost'+cls+'))')
status=execute('nprecls=float(n_elements(obsstrfe.strprepostfe'+cls+'.tdiff))')
status=execute('npreclshist=float(total(obsstrfe.yhistfepre'+cls+'))')
status=execute('yhistfeavgcls=(obsstrfe.yhistfepost'+cls+'+obsstrfe.yhistfepre'+cls+')/2.')

status=execute('nflpostcls=float(n_elements(obsstrfl.strprepostfl'+cls+'.ftdiff))')
status=execute('nflpostclshist=float(total(obsstrfl.yhistflpost'+cls+'))')
status=execute('nflprecls=float(n_elements(obsstrfl.strprepostfl'+cls+'.tdiff))')
status=execute('nflpreclshist=float(total(obsstrfl.yhistflpre'+cls+'))')
status=execute('yhistflavgcls=(obsstrfl.yhistflpost'+cls+'+obsstrfl.yhistflpre'+cls+')/2.')

;print,'obs FE is Ref. FL after FE AVG WT='+strtrim(total(obsstrfe.xhistfepost*obsstrfe.yhistfepostm5)/total(obsstrfe.yhistfepostm5),2)
;print,'obs FE is Ref. FL before FE AVG WT='+strtrim(total(obsstrfe.xhistfepre*obsstrfe.yhistfeprem5)/total(obsstrfe.yhistfeprem5),2)
;print,'obs FL is ref. FE after FL AVG WT='+strtrim(total(obsstrfl.xhistflpost*obsstrfl.yhistflpostm5)/total(obsstrfl.yhistflpostm5),2)
;print,'obs FL is ref. FE before FL AVG WT='+strtrim(total(obsstrfl.xhistflpre*obsstrfl.yhistflprem5)/total(obsstrfl.yhistflprem5),2)

status=execute('yhistfeprecls=obsstrfe.yhistfepre'+cls)
status=execute('yhistfepostcls=obsstrfe.yhistfepost'+cls)
status=execute('yhistflpostcls=obsstrfl.yhistflpost'+cls)
status=execute('yhistflprecls=obsstrfl.yhistflpre'+cls)

;yrobsfe=[0,1], yrobsfl=[0,0.1], yrsimfe=[0,1], yrsimfl=[0,0.1], yrboofe=[0,1], yrboofl=[0,0.1]}, $
;yrobsfed=[-1,1], yrobsfld=[-1,1], yrsimfed=[-1,1], yrsimfld=[-1,1], yrboofed=[-1,1], yrboofld=[-1,1]}


;Observed FE is reference for M5 flare
!p.multi=[12,2,6]
plot,[-1.5,obsstrfe.xhistfepost+1.5],[0,yhistfepostcls/npostcls],ps=10,yran=param.yrobsfe,ytit='N!dBIN!n / N!dTOT!n',ymarg=[0,3],xran=param.xran,xtit='',xtickname=strarr(10)+' ',ytickname=[' '],title='FE is Reference Event',/xsty,/ysty
sharpcorners
oploterr,obsstrfe.xhistfepost+1.5,yhistfepostcls/npostcls,sqrt(yhistfepostcls)/npostcls,ps=4
legend,['Matched Event After','Matched Event Before'],lines=[0,0],color=[0,180],/top,/right,chars=1
oplot,[-1.5,obsstrfe.xhistfepre+1.5],[0,yhistfeprecls/nprecls],ps=10,color=180
!p.color=180
oploterr,obsstrfe.xhistfepre+1.5,yhistfeprecls/nprecls,sqrt(yhistfeprecls)/nprecls,ps=4
!p.color=0
!p.multi=[10,2,6]
plot,[-1.5,obsstrfe.xhistfepost+1.5],[0,(yhistfepostcls-yhistfeprecls)/yhistfeavgcls],ytit='(N!dFE,'+cls+'!n - N!d'+cls+',FE!n) / Avg(N!d'+cls+',FE!n,N!dFE,'+cls+'!n)',ps=10,yran=param.yrobsfed,/ysty,xtit=textoidl('\tau')+'!dFE,'+cls+'!n [hrs]',xran=param.xran,ymarg=[3,0],/xsty
sharpcorners
hline,0,color=150
oploterr, obsstrfe.xhistfepost+1.5, (yhistfepostcls-yhistfeprecls)/yhistfeavgcls, fltarr(n_elements(obsstrfe.xhistfepost)), (sqrt(yhistfepostcls)+sqrt(yhistfeprecls))/yhistfeavgcls,ps=4,ymarg=2

xyouts,0.9,0.94,cls,/norm,color=0
xyouts,0.6,0.705,'Observed',/norm,color=0,chars=1.4
xyouts,0.6,0.38,'Poisson Sim.',/norm,color=0,chars=1.4
xyouts,0.6,0.055,'Re-ordered WTs (BS)',/norm,color=0,chars=1.4

xyouts,0.6,0.955,'A2',/norm,color=0,chars=1.4
xyouts,0.6,0.63,'B2',/norm,color=0,chars=1.4
xyouts,0.6,0.305,'C2',/norm,color=0,chars=1.4

xyouts,0.1,0.955,'A1',/norm,color=0,chars=1.4
xyouts,0.1,0.63,'B1',/norm,color=0,chars=1.4
xyouts,0.1,0.305,'C1',/norm,color=0,chars=1.4

;Observed Flare is reference
!p.multi=[11,2,6]
plot,[-1.5,obsstrfl.xhistflpost+1.5],[0,yhistflpostcls/nflpostcls],ps=10,yran=param.yrobsfl,ytit='N!dBIN!n / N!dTOT!n',ymarg=[0,3],xran=param.xran,xtit='',xtickname=strarr(10)+' ',ytickname=[' '],title='Flare is Reference Event',/xsty,/ysty
sharpcorners
oploterr,obsstrfl.xhistflpost+1.5,yhistflpostcls/nflpostcls,sqrt(yhistflpostcls)/nflpostcls,ps=4
;legend,['Flare Before','Flare After'],lines=[0,0],color=[0,180],/bottom,/right
oplot,[-1.5,obsstrfl.xhistflpre+1.5],[0,yhistflprecls/nflprecls],ps=10,color=180
!p.color=180
oploterr,obsstrfl.xhistflpre+1.5,yhistflprecls/nflprecls,sqrt(yhistflprecls)/nflprecls,ps=4
!p.color=0
!p.multi=[9,2,6]
plot,[-1.5,obsstrfl.xhistflpost+1.5],[0,(yhistflprecls-yhistflpostcls)/yhistflavgcls],ytit='(N!d'+cls+',FE!n - N!dFE,'+cls+'!n) / Avg(N!dFE,'+cls+'!n,N!d'+cls+',FE!n)',ps=10,yran=param.yrobsfld,/ysty,xtit=textoidl('\tau')+'!d'+cls+',FE!n [hrs]',xran=param.xran,ymarg=[3,0],/xsty
sharpcorners
hline,0,color=150
oploterr, obsstrfl.xhistflpost+1.5, (yhistflprecls-yhistflpostcls)/yhistflavgcls, fltarr(n_elements(obsstrfe.xhistfepost)), (sqrt(yhistfepostcls)+sqrt(yhistfeprecls))/yhistfeavgcls,ps=4,ymarg=2

;------------------------------------------------------------------->


status=execute('npostcls=float(n_elements(simstrfe.strprepostfe'+cls+'.ftdiff))')
status=execute('npostclshist=float(total(simstrfe.yhistfepost'+cls+'))')
status=execute('nprecls=float(n_elements(simstrfe.strprepostfe'+cls+'.tdiff))')
status=execute('npreclshist=float(total(simstrfe.yhistfepre'+cls+'))')
status=execute('yhistfeavgcls=(simstrfe.yhistfepost'+cls+'+simstrfe.yhistfepre'+cls+')/2.')

status=execute('nflpostcls=float(n_elements(simstrfl.strprepostfl'+cls+'.ftdiff))')
status=execute('nflpostclshist=float(total(simstrfl.yhistflpost'+cls+'))')
status=execute('nflprecls=float(n_elements(simstrfl.strprepostfl'+cls+'.tdiff))')
status=execute('nflpreclshist=float(total(simstrfl.yhistflpre'+cls+'))')
status=execute('yhistflavgcls=(simstrfl.yhistflpost'+cls+'+simstrfl.yhistflpre'+cls+')/2.')

;print,'sim FE is Ref. FL after FE AVG WT='+strtrim(total(simstrfe.xhistfepost*simstrfe.yhistfepostm5)/total(simstrfe.yhistfepostm5),2)
;print,'sim FE is Ref. FL before FE AVG WT='+strtrim(total(simstrfe.xhistfepre*simstrfe.yhistfeprem5)/total(simstrfe.yhistfeprem5),2)
;print,'sim FL is ref. FE after FL AVG WT='+strtrim(total(simstrfl.xhistflpost*simstrfl.yhistflpostm5)/total(simstrfl.yhistflpostm5),2)
;print,'sim FL is ref. FE before FL AVG WT='+strtrim(total(simstrfl.xhistflpre*simstrfl.yhistflprem5)/total(simstrfl.yhistflprem5),2)

status=execute('yhistfeprecls=simstrfe.yhistfepre'+cls)
status=execute('yhistfepostcls=simstrfe.yhistfepost'+cls)
status=execute('yhistflpostcls=simstrfl.yhistflpost'+cls)
status=execute('yhistflprecls=simstrfl.yhistflpre'+cls)

;Simulated FE is reference
!p.multi=[8,2,6]
plot,[-1.5,simstrfe.xhistfepost+1.5],[0,yhistfepostcls/npostcls],ps=10,yran=param.yrsimfe,ytit='N!dBIN!n / N!dTOT!n',ymarg=[1,2],xran=param.xran,xtit='',xtickname=strarr(10)+' ',ytickname=[' '],/xsty,/ysty
sharpcorners
oploterr,simstrfe.xhistfepost+1.5,yhistfepostcls/npostcls,sqrt(yhistfepostcls)/npostcls,ps=4
;legend,['Flare Before','Flare After'],lines=[0,0],color=[0,180],/bottom,/right
oplot,[-1.5,simstrfe.xhistfepre+1.5],[0,yhistfeprecls/nprecls],ps=10,color=180
!p.color=180
oploterr,simstrfe.xhistfepre+1.5,yhistfeprecls/nprecls,sqrt(yhistfeprecls)/nprecls,ps=4
!p.color=0
!p.multi=[6,2,6]
plot,[-1.5,simstrfe.xhistfepost+1.5],[0,(yhistfepostcls-yhistfeprecls)/yhistfeavgcls],ytit='(N!dFE,'+cls+'!n - N!d'+cls+',FE!n) / Avg(N!d'+cls+',FE!n,N!dFE,'+cls+'!n)',ps=10,yran=param.yrsimfed,/ysty,xtit=textoidl('\tau')+'!dFE,'+cls+'!n [hrs]',xran=param.xran,ymarg=[4,-1],/xsty
sharpcorners
hline,0,color=150
oploterr, simstrfe.xhistfepost+1.5, (yhistfepostcls-yhistfeprecls)/yhistfeavgcls, fltarr(n_elements(simstrfe.xhistfepost)), (sqrt(yhistfepostcls)+sqrt(yhistfeprecls))/yhistfeavgcls,ps=4,ymarg=2

;Simulated Flare is reference
!p.multi=[7,2,6]
plot,[-1.5,simstrfl.xhistflpost+1.5],[0,yhistflpostcls/nflpostcls],ps=10,yran=param.yrsimfl,ytit='N!dBIN!n / N!dTOT!n',ymarg=[1,2],xran=param.xran,xtit='',xtickname=strarr(10)+' ',ytickname=[' '],/xsty,/ysty
sharpcorners
oploterr,simstrfl.xhistflpost+1.5,yhistflpostcls/nflpostcls,sqrt(yhistflpostcls)/nflpostcls,ps=4
;legend,['Flare Before','Flare After'],lines=[0,0],color=[0,180],/bottom,/right
oplot,[-1.5,simstrfl.xhistflpre+1.5],[0,yhistflprecls/nflprecls],ps=10,color=180
!p.color=180
oploterr,simstrfl.xhistflpre+1.5,yhistflprecls/nflprecls,sqrt(yhistflprecls)/nflprecls,ps=4
!p.color=0
!p.multi=[5,2,6]
plot,[-1.5,simstrfl.xhistflpost+1.5],[0,(yhistflprecls-yhistflpostcls)/yhistflavgcls],ytit='(N!d'+cls+',FE!n - N!dFE,'+cls+'!n) / Avg(N!dFE,'+cls+'!n,N!d'+cls+',FE!n)',ps=10,yran=param.yrsimfld,/ysty,xtit=textoidl('\tau')+'!d'+cls+',FE!n [hrs]',xran=param.xran,ymarg=[4,-1],/xsty
sharpcorners
hline,0,color=150
oploterr, simstrfl.xhistflpost+1.5, (yhistflprecls-yhistflpostcls)/yhistflavgcls, fltarr(n_elements(simstrfe.xhistfepost)), (sqrt(yhistfepostcls)+sqrt(yhistfeprecls))/yhistfeavgcls,ps=4,ymarg=2



;------------------------------------------------------------------->

;Convert array of monte carlos to histograms
status=execute('yhistfepostcls=total(bootstr.yhistfepost'+cls+'arr,2)/param.nran')
status=execute('yhistfeprecls=total(bootstr.yhistfepre'+cls+'arr,2)/param.nran')

status=execute('yhistflpostcls=total(bootstr.yhistflpost'+cls+'arr,2)/param.nran')
status=execute('yhistflprecls=total(bootstr.yhistflpre'+cls+'arr,2)/param.nran')

;npostm5=float(n_elements(bootstr.strprepostfem5.ftdiff))
npostcls=float(total(yhistfepostcls))
;nprem5=float(n_elements(bootstr.strprepostfem5.tdiff))
nprecls=float(total(yhistfeprecls))
yhistfeavgcls=(yhistfepostcls+yhistfeprecls)/2.

;nflpostm5=float(n_elements(bootstr.strprepostflm5.ftdiff))
nflpostcls=float(total(yhistflpostcls))
;nflprem5=float(n_elements(bootstr.strprepostflm5.tdiff))
nflprecls=float(total(yhistflprecls))
yhistflavgcls=(yhistflpostcls+yhistflprecls)/2.

xhistfe=bootstr.xhistfe
xhistfl=bootstr.xhistfl

;print,'boot FE is Ref. FL after FE AVG WT='+strtrim(total(xhistfe*yhistfepostm5)/total(yhistfepostm5),2)
;print,'boot FE is Ref. FL before FE AVG WT='+strtrim(total(xhistfe*yhistfeprem5)/total(yhistfeprem5),2)
;print,'boot FL is ref. FE after FL AVG WT='+strtrim(total(xhistfl*yhistflpostm5)/total(yhistflpostm5),2)
;print,'boot FL is ref. FE before FL AVG WT='+strtrim(total(xhistfl*yhistflprem5)/total(yhistflprem5),2)



;Bootstrap FE is reference
;   YHISTFEPREM5ARR LONG      Array[57, 1000]
;   YHISTFEPREM1ARR LONG      Array[57, 1000]
;   YHISTFEPREC5ARR LONG      Array[57, 1000]
;   YHISTFEPREC1ARR LONG      Array[57, 1000]
;   YHISTFEPOSTM5ARR
;   YHISTFEPOSTM1ARR
;   YHISTFEPOSTC5ARR
;   YHISTFEPOSTC1ARR
;   YHISTFLPREM5ARR LONG      Array[57, 1000]
;   YHISTFLPREM1ARR LONG      Array[57, 1000]
;   YHISTFLPREC5ARR LONG      Array[57, 1000]
;   YHISTFLPREC1ARR LONG      Array[57, 1000]
;   YHISTFLPOSTM5ARR
;   YHISTFLPOSTM1ARR
;   YHISTFLPOSTC5ARR                                          
;   YHISTFLPOSTC1ARR                                          
;   XHISTFL         FLOAT     Array[57]                       
;   XHISTFE         FLOAT     Array[57]  

!p.multi=[4,2,6]
plot,[-1.5,xhistfe+1.5],[0,yhistfepostcls/npostcls],ps=10,yran=param.yrboofe,ytit='N!dBIN!n / N!dTOT!n',ymarg=[2,1],xran=param.xran,xtit='',xtickname=strarr(10)+' ',ytickname=[' '],/xsty,/ysty
sharpcorners
oploterr,xhistfe+1.5,yhistfepostcls/npostcls,sqrt(yhistfepostcls)/npostcls,ps=4
;legend,['Flare Before','Flare After'],lines=[0,0],color=[0,180],/bottom,/right
oplot,[-1.5,xhistfe+1.5],[0,yhistfeprecls/nprecls],ps=10,color=180
!p.color=180
oploterr,xhistfe+1.5,yhistfeprecls/nprecls,sqrt(yhistfeprecls)/nprecls,ps=4
!p.color=0
!p.multi=[2,2,6]
plot,[-1.5,xhistfe+1.5],[0,(yhistfepostcls-yhistfeprecls)/yhistfeavgcls],ytit='(N!dFE,'+cls+'!n - N!d'+cls+',FE!n) / Avg(N!d'+cls+',FE!n,N!dFE,'+cls+'!n)',ps=10,yran=param.yrboofed,/ysty,xtit=textoidl('\tau')+'!dFE,'+cls+'!n [hrs]',xran=param.xran,ymarg=[5,-2],/xsty
sharpcorners
hline,0,color=150
oploterr, xhistfe+1.5, (yhistfepostcls-yhistfeprecls)/yhistfeavgcls, fltarr(n_elements(xhistfe)), (sqrt(yhistfepostcls)+sqrt(yhistfeprecls))/yhistfeavgcls,ps=4,ymarg=2

;Bootstrap Flare is reference
!p.multi=[3,2,6]
plot,[-1.5,xhistfl+1.5],[0,yhistflpostcls/nflpostcls],ps=10,yran=param.yrboofl,ytit='N!dBIN!n / N!dTOT!n',ymarg=[2,1],xran=param.xran,xtit='',xtickname=strarr(10)+' ',ytickname=[' '],/xsty,/ysty
sharpcorners
oploterr,xhistfl+1.5,yhistflpostcls/nflpostcls,sqrt(yhistflpostcls)/nflpostcls,ps=4
;legend,['Flare Before','Flare After'],lines=[0,0],color=[0,180],/bottom,/right
oplot,[-1.5,xhistfl+1.5],[0,yhistflprecls/nflprecls],ps=10,color=180
!p.color=180
oploterr,xhistfl+1.5,yhistflprecls/nflprecls,sqrt(yhistflprecls)/nflprecls,ps=4
!p.color=0
!p.multi=[1,2,6]
plot,[-1.5,xhistfl+1.5],[0,(yhistflprecls-yhistflpostcls)/yhistflavgcls],ytit='(N!d'+cls+',FE!n - N!dFE,'+cls+'!n) / Avg(N!dFE,'+cls+'!n,N!d'+cls+',FE!n)',ps=10,yran=param.yrboofld,/ysty,xtit=textoidl('\tau')+'!d'+cls+',FE!n [hrs]',xran=param.xran,ymarg=[5,-2],/xsty
sharpcorners
hline,0,color=150
oploterr, xhistfl+1.5, (yhistflprecls-yhistflpostcls)/yhistflavgcls, fltarr(n_elements(xhistfe)), (sqrt(yhistfepostcls)+sqrt(yhistfeprecls))/yhistfeavgcls,ps=4,ymarg=2


closeplotenv

stop

end

;----------------------------------------------------------------------------->

;Loop over the 4 flare thresholds

pro fefl_sixpanel_generalize

root='~/science/projects/stereo_sympathetic_flaring/'
datapath=root+'data/'
plotpath=root+'plots/determine_fe-fl_pair_relation/png/'
epspath=root+'plots/determine_fe-fl_pair_relation/eps/'

;Restore observed lists
restore,datapath+'determine_fe-fl_pair_relation.res1.sav',/ver
restore,datapath+'determine_fe-fl_pair_relation_festruct.sav',/ver
restore,datapath+'determine_fe-fl_pair_relation_flstruct.sav',/ver

;Observed timeseries
tlec1 = timle[wle[wc1]]
tlec5 = timle[wle[wc5]]
tleM1 = timle[wle[wM1]]
tleM5 = timle[wle[wM5]]
tfe = timfe[wfe[wfe2]]

;Restore simulated lists
restore,datapath+'simulate_fe_fl_pair_relation-festruct.sav',/ver
restore,datapath+'simulate_fe_fl_pair_relation-flstruct.sav',/ver

;Restore bootstrapped lists
restore,datapath+'simulate_fe_fl_pair_bootstrap-struct.sav',/ver

stop

param={cls:'M5',xran:[-2.,24.],epspath:epspath,nran:1000., $
yrobsfe:[0,0.02], yrobsfl:[0,0.2], yrsimfe:[0,0.02], yrsimfl:[0,0.2], yrboofe:[0,0.06], yrboofl:[0,0.15], $
yrobsfed:[-2,2], yrobsfld:[-3,3], yrsimfed:[-2,2], yrsimfld:[-3,3], yrboofed:[-2,2], yrboofld:[-2,2]}
fefl_sixpanel_plotdist,param,obsstrfe,obsstrfl,simstrfe,simstrfl,bootstr

param={cls:'M1',xran:[-2.,24.],epspath:epspath,nran:1000., $ 
yrobsfe:[0,0.05], yrobsfl:[0,0.1], yrsimfe:[0,0.05], yrsimfl:[0,0.15], yrboofe:[0,0.07], yrboofl:[0,0.1], $
yrobsfed:[-1,1], yrobsfld:[-1,1], yrsimfed:[-1,1], yrsimfld:[-1,1], yrboofed:[-1,1], yrboofld:[-1,1]}
fefl_sixpanel_plotdist,param,obsstrfe,obsstrfl,simstrfe,simstrfl,bootstr

param={cls:'C5',xran:[-2.,24.],epspath:epspath,nran:1000., $ 
yrobsfe:[0,1], yrobsfl:[0,0.1], yrsimfe:[0,1], yrsimfl:[0,0.15], yrboofe:[0,1], yrboofl:[0,0.1], $
yrobsfed:[-1,1], yrobsfld:[-1,1], yrsimfed:[-1,1], yrsimfld:[-1,1], yrboofed:[-1,1], yrboofld:[-1,1]}
fefl_sixpanel_plotdist,param,obsstrfe,obsstrfl,simstrfe,simstrfl,bootstr

param={cls:'C1',xran:[-2.,24.],epspath:epspath,nran:1000., $
yrobsfe:[0,1], yrobsfl:[0,0.1], yrsimfe:[0,1], yrsimfl:[0,0.1], yrboofe:[0,1], yrboofl:[0,0.1], $
yrobsfed:[-1,1], yrobsfld:[-1,1], yrsimfed:[-1,1], yrsimfld:[-1,1], yrboofed:[-1,1], yrboofld:[-1,1]}
fefl_sixpanel_plotdist,param,obsstrfe,obsstrfl,simstrfe,simstrfl,bootstr



stop

;Determine the KS probs for each pair of distributions
ksobsstr={dfem5:0.,pfem5:0.,dfem1:0.,pfem1:0.,dfec5:0.,pfec5:0.,dfec1:0.,pfec1:0., $
	dflm5:0.,pflm5:0.,dflm1:0.,pflm1:0.,dflc5:0.,pflc5:0.,dflc1:0.,pflc1:0.}

ksobsstr48={dfem5:0.,pfem5:0.,dfem1:0.,pfem1:0.,dfec5:0.,pfec5:0.,dfec1:0.,pfec1:0., $
	dflm5:0.,pflm5:0.,dflm1:0.,pflm1:0.,dflc5:0.,pflc5:0.,dflc1:0.,pflc1:0.}

ksbootstr={dfem5:0.,pfem5:0.,dfem1:0.,pfem1:0.,dfec5:0.,pfec5:0.,dfec1:0.,pfec1:0., $
	dflm5:0.,pflm5:0.,dflm1:0.,pflm1:0.,dflc5:0.,pflc5:0.,dflc1:0.,pflc1:0.}

ksbootstr48={dfem5:0.,pfem5:0.,dfem1:0.,pfem1:0.,dfec5:0.,pfec5:0.,dfec1:0.,pfec1:0., $
	dflm5:0.,pflm5:0.,dflm1:0.,pflm1:0.,dflc5:0.,pflc5:0.,dflc1:0.,pflc1:0.}

kssimstr={dfem5:0.,pfem5:0.,dfem1:0.,pfem1:0.,dfec5:0.,pfec5:0.,dfec1:0.,pfec1:0., $
	dflm5:0.,pflm5:0.,dflm1:0.,pflm1:0.,dflc5:0.,pflc5:0.,dflc1:0.,pflc1:0.}

kssimstr48={dfem5:0.,pfem5:0.,dfem1:0.,pfem1:0.,dfec5:0.,pfec5:0.,dfec1:0.,pfec1:0., $
	dflm5:0.,pflm5:0.,dflm1:0.,pflm1:0.,dflc5:0.,pflc5:0.,dflc1:0.,pflc1:0.}

stop

;Observations 0-24 hours
kstwo, (obsstrfe.yhistfeprem5)[0:8], (obsstrfe.yhistfepostm5)[0:8], dfem5, pfem5
ksobsstr.dfem5=dfem5 & ksobsstr.pfem5=pfem5
kstwo, (obsstrfe.yhistfeprem1)[0:8], (obsstrfe.yhistfepostm1)[0:8], dfem1, pfem1
ksobsstr.dfem1=dfem1 & ksobsstr.pfem1=pfem1
kstwo, (obsstrfe.yhistfeprec5)[0:8], (obsstrfe.yhistfepostc5)[0:8], dfec5, pfec5
ksobsstr.dfec5=dfec5 & ksobsstr.pfec5=pfec5
kstwo, (obsstrfe.yhistfeprec1)[0:8], (obsstrfe.yhistfepostc1)[0:8], dfec1, pfec1
ksobsstr.dfec1=dfec1 & ksobsstr.pfec1=pfec1

kstwo, (obsstrfl.yhistflprem5)[0:8], (obsstrfl.yhistflpostm5)[0:8], dflm5, pflm5
ksobsstr.dflm5=dflm5 & ksobsstr.pflm5=pflm5
kstwo, (obsstrfl.yhistflprem1)[0:8], (obsstrfl.yhistflpostm1)[0:8], dflm1, pflm1
ksobsstr.dflm1=dflm1 & ksobsstr.pflm1=pflm1
kstwo, (obsstrfl.yhistflprec5)[0:8], (obsstrfl.yhistflpostc5)[0:8], dflc5, pflc5
ksobsstr.dflc5=dflc5 & ksobsstr.pflc5=pflc5
kstwo, (obsstrfl.yhistflprec1)[0:8], (obsstrfl.yhistflpostc1)[0:8], dflc1, pflc1
ksobsstr.dflc1=dflc1 & ksobsstr.pflc1=pflc1

;Observations 0-48 hours
kstwo, (obsstrfe.yhistfeprem5)[0:16], (obsstrfe.yhistfepostm5)[0:16], dfem5, pfem5
ksobsstr48.dfem5=dfem5 & ksobsstr48.pfem5=pfem5
kstwo, (obsstrfe.yhistfeprem1)[0:16], (obsstrfe.yhistfepostm1)[0:16], dfem1, pfem1
ksobsstr48.dfem1=dfem1 & ksobsstr48.pfem1=pfem1
kstwo, (obsstrfe.yhistfeprec5)[0:16], (obsstrfe.yhistfepostc5)[0:16], dfec5, pfec5
ksobsstr48.dfec5=dfec5 & ksobsstr48.pfec5=pfec5
kstwo, (obsstrfe.yhistfeprec1)[0:16], (obsstrfe.yhistfepostc1)[0:16], dfec1, pfec1
ksobsstr48.dfec1=dfec1 & ksobsstr48.pfec1=pfec1

kstwo, (obsstrfl.yhistflprem5)[0:16], (obsstrfl.yhistflpostm5)[0:16], dflm5, pflm5
ksobsstr48.dflm5=dflm5 & ksobsstr48.pflm5=pflm5
kstwo, (obsstrfl.yhistflprem1)[0:16], (obsstrfl.yhistflpostm1)[0:16], dflm1, pflm1
ksobsstr48.dflm1=dflm1 & ksobsstr48.pflm1=pflm1
kstwo, (obsstrfl.yhistflprec5)[0:16], (obsstrfl.yhistflpostc5)[0:16], dflc5, pflc5
ksobsstr48.dflc5=dflc5 & ksobsstr48.pflc5=pflc5
kstwo, (obsstrfl.yhistflprec1)[0:16], (obsstrfl.yhistflpostc1)[0:16], dflc1, pflc1
ksobsstr48.dflc1=dflc1 & ksobsstr48.pflc1=pflc1

;Simulations 0-24 hours
kstwo, (simstrfe.yhistfeprem5)[0:8], (simstrfe.yhistfepostm5)[0:8], dfem5, pfem5
kssimstr.dfem5=dfem5 & kssimstr.pfem5=pfem5
kstwo, (simstrfe.yhistfeprem1)[0:8], (simstrfe.yhistfepostm1)[0:8], dfem1, pfem1
kssimstr.dfem1=dfem1 & kssimstr.pfem1=pfem1
kstwo, (simstrfe.yhistfeprec5)[0:8], (simstrfe.yhistfepostc5)[0:8], dfec5, pfec5
kssimstr.dfec5=dfec5 & kssimstr.pfec5=pfec5
kstwo, (simstrfe.yhistfeprec1)[0:8], (simstrfe.yhistfepostc1)[0:8], dfec1, pfec1
kssimstr.dfec1=dfec1 & kssimstr.pfec1=pfec1

kstwo, (simstrfl.yhistflprem5)[0:8], (simstrfl.yhistflpostm5)[0:8], dflm5, pflm5
kssimstr.dflm5=dflm5 & kssimstr.pflm5=pflm5
kstwo, (simstrfl.yhistflprem1)[0:8], (simstrfl.yhistflpostm1)[0:8], dflm1, pflm1
kssimstr.dflm1=dflm1 & kssimstr.pflm1=pflm1
kstwo, (simstrfl.yhistflprec5)[0:8], (simstrfl.yhistflpostc5)[0:8], dflc5, pflc5
kssimstr.dflc5=dflc5 & kssimstr.pflc5=pflc5
kstwo, (simstrfl.yhistflprec1)[0:8], (simstrfl.yhistflpostc1)[0:8], dflc1, pflc1
kssimstr.dflc1=dflc1 & kssimstr.pflc1=pflc1

;Simulations 0-48 hours
kstwo, (simstrfe.yhistfeprem5)[0:16], (simstrfe.yhistfepostm5)[0:16], dfem5, pfem5
kssimstr48.dfem5=dfem5 & kssimstr48.pfem5=pfem5
kstwo, (simstrfe.yhistfeprem1)[0:16], (simstrfe.yhistfepostm1)[0:16], dfem1, pfem1
kssimstr48.dfem1=dfem1 & kssimstr48.pfem1=pfem1
kstwo, (simstrfe.yhistfeprec5)[0:16], (simstrfe.yhistfepostc5)[0:16], dfec5, pfec5
kssimstr48.dfec5=dfec5 & kssimstr48.pfec5=pfec5
kstwo, (simstrfe.yhistfeprec1)[0:16], (simstrfe.yhistfepostc1)[0:16], dfec1, pfec1
kssimstr48.dfec1=dfec1 & kssimstr48.pfec1=pfec1

kstwo, (simstrfl.yhistflprem5)[0:16], (simstrfl.yhistflpostm5)[0:16], dflm5, pflm5
kssimstr48.dflm5=dflm5 & kssimstr48.pflm5=pflm5
kstwo, (simstrfl.yhistflprem1)[0:16], (simstrfl.yhistflpostm1)[0:16], dflm1, pflm1
kssimstr48.dflm1=dflm1 & kssimstr48.pflm1=pflm1
kstwo, (simstrfl.yhistflprec5)[0:16], (simstrfl.yhistflpostc5)[0:16], dflc5, pflc5
kssimstr48.dflc5=dflc5 & kssimstr48.pflc5=pflc5
kstwo, (simstrfl.yhistflprec1)[0:16], (simstrfl.yhistflpostc1)[0:16], dflc1, pflc1
kssimstr48.dflc1=dflc1 & kssimstr48.pflc1=pflc1



;Need to loop over 1000 hists and calc. KS for each pair
fefl_sixpanel_generalize_bootks,bootstr,ksbootstr,xran=[0,8]

fefl_sixpanel_generalize_bootks,bootstr,ksbootstr48,xran=[0,16]


;Bootstrapped 0-24 hours
;kstwo, (bootstr.yhistfeprem5arr)[0:8], (bootstr.yhistfepostm5arr)[0:8], dfem5, pfem5
;ksbootstr.dfem5=dfem5 & ksbootstr.pfem5=pfem5
;kstwo, (bootstr.yhistfeprem1arr)[0:8], (bootstr.yhistfepostm1arr)[0:8], dfem1, pfem1
;ksbootstr.dfem1=dfem1 & ksbootstr.pfem1=pfem1
;kstwo, (bootstr.yhistfeprec5arr)[0:8], (bootstr.yhistfepostc5arr)[0:8], dfec5, pfec5
;ksbootstr.dfec5=dfec5 & ksbootstr.pfec5=pfec5
;kstwo, (bootstr.yhistfeprec1arr)[0:8], (bootstr.yhistfepostc1arr)[0:8], dfec1, pfec1
;ksbootstr.dfec1=dfec1 & ksbootstr.pfec1=pfec1

;kstwo, (bootstr.yhistflprem5arr)[0:8], (bootstr.yhistflpostm5arr)[0:8], dflm5, pflm5
;ksbootstr.dflm5=dflm5 & ksbootstr.pflm5=pflm5
;kstwo, (bootstr.yhistflprem1arr)[0:8], (bootstr.yhistflpostm1arr)[0:8], dflm1, pflm1
;ksbootstr.dflm1=dflm1 & ksbootstr.pflm1=pflm1
;kstwo, (bootstr.yhistflprec5arr)[0:8], (bootstr.yhistflpostc5arr)[0:8], dflc5, pflc5
;ksbootstr.dflc5=dflc5 & ksbootstr.pflc5=pflc5
;kstwo, (bootstr.yhistflprec1arr)[0:8], (bootstr.yhistflpostc1arr)[0:8], dflc1, pflc1
;ksbootstr.dflc1=dflc1 & ksbootstr.pflc1=pflc1

;Bootstrapped 0-48 hours
;kstwo, (bootstr.yhistfeprem5arr)[0:16], (bootstr.yhistfepostm5arr)[0:16], dfem5, pfem5
;ksbootstr48.dfem5=dfem5 & ksbootstr48.pfem5=pfem5
;kstwo, (bootstr.yhistfeprem1arr)[0:16], (bootstr.yhistfepostm1arr)[0:16], dfem1, pfem1
;ksbootstr48.dfem1=dfem1 & ksbootstr48.pfem1=pfem1
;kstwo, (bootstr.yhistfeprec5arr)[0:16], (bootstr.yhistfepostc5arr)[0:16], dfec5, pfec5
;ksbootstr48.dfec5=dfec5 & ksbootstr48.pfec5=pfec5
;kstwo, (bootstr.yhistfeprec1arr)[0:16], (bootstr.yhistfepostc1arr)[0:16], dfec1, pfec1
;ksbootstr48.dfec1=dfec1 & ksbootstr48.pfec1=pfec1

;kstwo, (bootstr.yhistflprem5arr)[0:16], (bootstr.yhistflpostm5arr)[0:16], dflm5, pflm5
;ksbootstr48.dflm5=dflm5 & ksbootstr48.pflm5=pflm5
;kstwo, (bootstr.yhistflprem1arr)[0:16], (bootstr.yhistflpostm1arr)[0:16], dflm1, pflm1
;ksbootstr48.dflm1=dflm1 & ksbootstr48.pflm1=pflm1
;kstwo, (bootstr.yhistflprec5arr)[0:16], (bootstr.yhistflpostc5arr)[0:16], dflc5, pflc5
;ksbootstr48.dflc5=dflc5 & ksbootstr48.pflc5=pflc5
;kstwo, (bootstr.yhistflprec1arr)[0:16], (bootstr.yhistflpostc1arr)[0:16], dflc1, pflc1
;ksbootstr48.dflc1=dflc1 & ksbootstr48.pflc1=pflc1

tablearr=strarr(5,8)
tablearr[1,*]=['FE/M5','FE/M1','FE/C5','FE/C1','M5/FE','M1/FE','C5/FE','C1/FE']
tablearr[2,*]=string([ksobsstr.pfem5,ksobsstr.pfem1,ksobsstr.pfec5,ksobsstr.pfec1,ksobsstr.pflm5,ksobsstr.pflm1,ksobsstr.pflc5,ksobsstr.pflc1],form='(F6.3)')
tablearr[3,*]=string([kssimstr.pfem5,kssimstr.pfem1,kssimstr.pfec5,kssimstr.pfec1,kssimstr.pflm5,kssimstr.pflm1,kssimstr.pflc5,kssimstr.pflc1],form='(F6.3)')
tablearr[4,*]=string([ksbootstr.pfem5,ksbootstr.pfem1,ksbootstr.pfec5,ksbootstr.pfec1,ksbootstr.pflm5,ksbootstr.pflm1,ksbootstr.pflc5,ksbootstr.pflc1],form='(F6.3)')
array2textable, tablearr, ['Fig. \#', 'Ref. Event', 'Obs.', 'Sim.', 'Bootstrap'], 'this is a table about blah', label='tableblah', outfile=root+'paper/ks_table24.txt'

tablearr48=strarr(5,8)
tablearr48[1,*]=['FE/M5','FE/M1','FE/C5','FE/C1','M5/FE','M1/FE','C5/FE','C1/FE']
tablearr48[2,*]=string([ksobsstr48.pfem5,ksobsstr48.pfem1,ksobsstr48.pfec5,ksobsstr48.pfec1,ksobsstr48.pflm5,ksobsstr48.pflm1,ksobsstr48.pflc5,ksobsstr48.pflc1],form='(F6.3)')
tablearr48[3,*]=string([kssimstr48.pfem5,kssimstr48.pfem1,kssimstr48.pfec5,kssimstr48.pfec1,kssimstr48.pflm5,kssimstr48.pflm1,kssimstr48.pflc5,kssimstr48.pflc1],form='(F6.3)')
tablearr48[4,*]=string([ksbootstr48.pfem5,ksbootstr48.pfem1,ksbootstr48.pfec5,ksbootstr48.pfec1,ksbootstr48.pflm5,ksbootstr48.pflm1,ksbootstr48.pflc5,ksbootstr48.pflc1],form='(F6.3)')
array2textable, tablearr48, ['Fig. \#', 'Ref. Event', 'Obs.', 'Sim.', 'Bootstrap'], 'this is a table for 48hr time range', label='tableblah48', outfile=root+'paper/ks_table48.txt'

stop

print,[ksobsstr.pfem5,ksobsstr.pfem1,ksobsstr.pfec5,ksobsstr.pfec1,ksobsstr.pflm5,ksobsstr.pflm1,ksobsstr.pflc5,ksobsstr.pflc1]

print,[ksobsstr48.pfem5,ksobsstr48.pfem1,ksobsstr48.pfec5,ksobsstr48.pfec1,ksobsstr48.pflm5,ksobsstr48.pflm1,ksobsstr48.pflc5,ksobsstr48.pflc1]

print,[kssimstr.pfem5,kssimstr.pfem1,kssimstr.pfec5,kssimstr.pfec1,kssimstr.pflm5,kssimstr.pflm1,kssimstr.pflc5,kssimstr.pflc1]

print,[kssimstr48.pfem5,kssimstr48.pfem1,kssimstr48.pfec5,kssimstr48.pfec1,kssimstr48.pflm5,kssimstr48.pflm1,kssimstr48.pflc5,kssimstr48.pflc1]

stop

restore,datapath+'determine_fe-fl_pair_relation_res1.sav',/ver

;Check the times of the flare-FE pairs in the first bin of the MF histograms 

wfetdiff=where(abs(strprepostfem5.tdiff) le 3.*3600.)

print,anytim((strprepostfem5.t1)[wfetdiff],/vms)
print,anytim((strprepostfem5.t2)[wfetdiff],/vms)

wfeftdiff=where(abs(strprepostfem5.ftdiff) le 3.*3600.)

print,anytim((strprepostfem5.ft1)[wfeftdiff],/vms) 
print,anytim((strprepostfem5.ft2)[wfeftdiff],/vms)

wfltdiff=where(abs(strprepostflm5.tdiff) le 3.*3600.)

print,anytim((strprepostflm5.t1)[wfltdiff],/vms)
print,anytim((strprepostflm5.t2)[wfltdiff],/vms)

wflftdiff=where(abs(strprepostflm5.ftdiff) le 3.*3600.)

print,anytim((strprepostflm5.ft1)[wflftdiff],/vms) 
print,anytim((strprepostflm5.ft2)[wflftdiff],/vms)

stop

end
