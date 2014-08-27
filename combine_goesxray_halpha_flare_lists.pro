;2014-04-17: combine the NGDC and GOES flare lists. There seems to be many large flares in GOES list with no positions, and many NGDC Halpha flares with positions, but no goes class.
;Attempting to combine the two to get a more complete flare list

;This code requires SSW and this repository: https://github.com/pohuigin/gen_library

pro combine_goesxray_halpha_flare_lists,res1tore=res1tore

dpath='~/science/data/helio_flare_lists/'
fgoes=dpath+['hec_goes_list_1996-1997.txt','hec_goes_list_1998-1999.txt','hec_goes_list_2000-2001.txt','hec_goes_list_2002-2003.txt','hec_goes_list_2004-2005.txt','hec_goes_list_2006-2007.txt','hec_goes_list_2008-2009.txt','hec_goes_list_2010-2011.txt','hec_goes_list_2012-2013.txt','hec_goes_list_2014.txt']
fhalpha=dpath+['hec_halpha_list_1996-1997.txt','hec_halpha_list_1998-1999.txt','hec_halpha_list_2000-2001.txt','hec_halpha_list_2002-2003.txt','hec_halpha_list_2004-2005.txt','hec_halpha_list_2006-2007.txt','hec_halpha_list_2008-2009.txt','hec_halpha_list_2010.txt']

savpath='~/science/projects/stereo_sympathetic_flaring/data/'

if keyword_set(res1tore) then goto, skiptores1

for i=0l,n_elements(fhalpha)-1l do begin
	readcol,fhalpha[i],halphaline,delim='$%',form='A'
	if i eq 0 then halphalinearr=halphaline else halphalinearr=[halphalinearr,halphaline]
endfor
halphaline=halphalinearr


for i=0l,n_elements(fgoes)-1l do begin
	readcol,fgoes[i],goesline,delim='$%',form='A'
	if i eq 0 then goeslinearr=goesline else goeslinearr=[goeslinearr,goesline]
endfor
goesline=goeslinearr



ngoes=n_elements(goesline)
nhalpha=n_elements(halphaline)

goesarr=strarr(9,ngoes)
halphaarr=strarr(9,nhalpha)

for i=0l,ngoes-1l do goesarr[*,i]=(str_sep(goesline[i],' '))[0:8]

for i=0l,nhalpha-1l do halphaarr[*,i]=(str_sep(halphaline[i],' '))[0:8]




;Check for valid dates in the Halpha list
hyyyy=strmid(halphaarr[1,1:*],0,4)
hwgood=where(hyyyy ge 1996 and hyyyy le 2014)

;Check for valid dates in the GOES list
gyyyy=strmid(goesarr[1,1:*],0,4)
gwgood=where(gyyyy ge 1996 and gyyyy le 2014)

htime=anytim((halphaarr[1,1:*])[hwgood],/vms)
gtime=anytim((goesarr[1,1:*])[gwgood],/vms)




;Update the arrays so there is one less element (leave out the header)
ngoes=n_elements(gtime)
nhalpha=n_elements(htime)

whpeak=where((halphaarr[2,1:*])[hwgood] ne '""')
wnohpeak=where((halphaarr[2,1:*])[hwgood] eq '""')

htim0=anytim((halphaarr[1,1:*])[hwgood])
htimp=lonarr(nhalpha)
htimp[whpeak]=anytim(((halphaarr[2,1:*])[hwgood])[whpeak])
htim1=anytim((halphaarr[3,1:*])[hwgood])

;Fill missing halpha peak times with the average of start and end
htimp[wnohpeak]=round((htim1[wnohpeak]+htim0[wnohpeak])/2.)

gtim0=anytim((goesarr[1,1:*])[gwgood])
gtimp=anytim((goesarr[2,1:*])[gwgood])
gtim1=anytim((goesarr[3,1:*])[gwgood])

hlon=(halphaarr[6,1:*])[hwgood]
whlon=where(hlon ne '""')
wnohlon=where(hlon eq '""')
hlon[wnohlon]=1000
hlon=float(hlon)

glon=(goesarr[6,1:*])[gwgood]
wglon=where(glon ne '""')
wnoglon=where(glon eq '""')
glon[wnoglon]=1000
glon=float(glon)

hlat=float((halphaarr[5,1:*])[hwgood])
hlat[wnohlon]=1000

glat=float((goesarr[5,1:*])[gwgood])
glat[wnoglon]=1000

hclass=(halphaarr[8,1:*])[hwgood]
gclass=(goesarr[8,1:*])[gwgood]

hnoaa=(halphaarr[4,1:*])[hwgood]
gnoaa=(goesarr[4,1:*])[gwgood]

;MAKE A FLARE LIST FOR ALAN
;save,gtimp,gclass,file=dpath+'alan_all_goes_tim_class.sav'
;stop
;gtimp=anytim((goesarr[2,1:*])[gwgood])
;gclass=(goesarr[8,1:*])[gwgood]
;walan=where(gtimp ge anytim('1-jun-1997') and gtimp le anytim('1-jun-2008'))
;alanstr={time:'',class:''}
;alanstr=replicate(alanstr,n_elements(walan))
;alanstr.time=anytim(gtimp[walan],/ecs)
;alanstr.class=gclass[walan]
;write_data_file,alanstr,file=dpath+'alan_all_goes_19970601_20080601.txt'



;loop over all goes events with missing locations and try to match each one with an halpha event
nnogoes=n_elements(wnoglon)

;Array of matched positions in the halpha array
match_hind=fltarr(nnogoes)-1.
match_gind=fltarr(nnogoes)-1.

;Use only the halpha event with positions for matching
htimpm=htimp[whlon]


for i=0l,nnogoes-1l do begin
	thiswg=wnoglon[i]
	
	wmatch=where(htimpm ge gtim0[thiswg] and htimpm le gtim1[thiswg])
	
	if n_elements(wmatch) gt 1 then wmatch=(where(abs(htimpm-gtimp[thiswg]) eq min(abs(htimpm-gtimp[thiswg]))))[0]
	
;	if wmatch[0] eq -1 then begin
;		wmatch=where(htimpm ge (gtim0[thiswg]-3600.*2.) and htimpm le (gtim1[thiswg]+3600.*2.))
;		if n_elements(wmatch) gt 1 then wmatch=(where(abs(htimpm-gtimp[thiswg]) eq min(abs(htimpm-gtimp[thiswg]))))[0]
;	endif
	
	match_hind[i]=wmatch
	match_gind[i]=thiswg
endfor

help,match_hind

;Find where matches were successful in the no-goes-pos. where array
wgoodmatch=where(match_hind ne -1)

;Only take good matches in the halpha and goes where no-goes-pos. arrays
wnoglon=wnoglon[wgoodmatch]
match_hind=match_hind[wgoodmatch]
match_gind=match_gind[wgoodmatch]

;Fill structure with combined info
flrcombstr={gst:0l,gpk:0l,gen:0l,hst:-1l,hpk:-1l,hen:-1l,gclass:'""',hclass:'""',gnoaa:-1,hnoaa:-1,lat:1000.,lon:1000.,halphamatch:-1,narmatch:-1}

flrcombstr=replicate(flrcombstr,ngoes)

;Goes info
flrcombstr.gst=reform(gtim0)
flrcombstr.gpk=reform(gtimp)
flrcombstr.gen=reform(gtim1)
flrcombstr.gclass=reform(gclass)
flrcombstr.gnoaa=reform(gnoaa)

flrcombstr.lat=reform(glat)
flrcombstr.lon=reform(glon)

flrcombstr[wglon].halphamatch=0

;Halpha info for matched events
flrcombstr[wnoglon].hst=htim0[whlon[match_hind]]
flrcombstr[wnoglon].hpk=htimp[whlon[match_hind]]
flrcombstr[wnoglon].hen=htim1[whlon[match_hind]]
flrcombstr[wnoglon].hclass=hclass[whlon[match_hind]]
flrcombstr[wnoglon].hnoaa=hnoaa[whlon[match_hind]]

flrcombstr[wnoglon].lat=hlat[whlon[match_hind]]
flrcombstr[wnoglon].lon=hlon[whlon[match_hind]]

flrcombstr[wnoglon].halphamatch=1


;Fill in the missing goes NOAA #s with matched Halpha
wgmissnoaa=where(flrcombstr.gnoaa le 0)
whnoaa=where(flrcombstr[wgmissnoaa].hnoaa gt 0)
flrcombstr[wgmissnoaa[whnoaa]].gnoaa=flrcombstr[wgmissnoaa[whnoaa]].hnoaa

save,flrcombstr,goesarr,halphaarr,halphaline,goesline,file=savpath+'combine_goesxray_halpha_flare_lists-res1.sav'




;Loop through and for every flare with no position, but has a noaa number, fill in the position based on the RD_NAR info for that day
whasnoaa=where(flrcombstr.gnoaa gt 0)
wnopos=where(flrcombstr[whasnoaa].lat eq 1000)


nnoaanopos=n_elements(wnopos)

for i=0,nnoaanopos-1 do begin

	flrcombstr[whasnoaa[wnopos[i]]].narmatch=0
	thisnoaa=flrcombstr[whasnoaa[wnopos[i]]].gnoaa

	thisgpk=flrcombstr[whasnoaa[wnopos[i]]].gpk

	thisdate0=anytim(flrcombstr[whasnoaa[wnopos[i]]].gpk-3600.*49.,/vms,/date)
	thisdate1=anytim(flrcombstr[whasnoaa[wnopos[i]]].gpk+3600.*49.,/vms,/date)
	rd_nar,thisdate0+' 00:00',thisdate1+' 23:59',narstr;,/nearest

if data_type(narstr) ne 8 then continue

	wnoaalt1000=where(narstr.noaa lt 5000)
	if (wnoaalt1000)[0] ne -1 then narstr[wnoaalt1000].noaa=narstr[wnoaalt1000].noaa+10000

	wnarstr=(where(narstr.noaa eq thisnoaa))[0]

;if strmid(flrcombstr[whasnoaa[wnopos[i]]].gclass,0,1) eq 'X' then stop

	if wnarstr[0] ne -1 then begin
	
		thisnar=narstr[wnarstr[0]]
		nartim=anytim(thisnar)
		
		thislonlat=thisnar.location

		thisnardlon=diff_rot((thisgpk-nartim)/3600./24.,thislonlat[1])
	
;Include the differential rotation btwn nar time and flr time
		flrcombstr[whasnoaa[wnopos[i]]].lon=thislonlat[0]+thisnardlon
		flrcombstr[whasnoaa[wnopos[i]]].lat=thislonlat[1]
		flrcombstr[whasnoaa[wnopos[i]]].narmatch=1
	endif else begin
	
if strmid(flrcombstr[whasnoaa[wnopos[i]]].gclass,0,1) eq 'X' then stop

;		rd_nar,thisdate+' 00:00',thisdate+' 23:59',gevstr,/nearest
		
		
		
		
	endelse

endfor





;WHY CAN't I MATCH THE X-CLASS FLARES?!?!

;FOR SOME REASON ITS NOT LOADING THEM IN BUT ALSO ISNT STOPPING THE CODE!!




;matchstr=match_flares(gtime[wglon], htime, glon[wglon], hlon, glat, hlat, gclass, hclass, $
;	timthresh=3600.,radthresh=10., debug=0, match12=match12, match21=match21, label1='goes', label2='halpha', doplot=1, bookkeeping=bookstr)








utplot,flrcombstr.gpk,flrcombstr.lat,ps=4,yran=[-100,100]
setcolors,/sys
oplot,flrcombstr[where(flrcombstr.HALPHAMATCH eq 1)].gpk,flrcombstr[where(flrcombstr.HALPHAMATCH eq 1)].lat,ps=4,color=!red







save,flrcombstr,goesarr,halphaarr,halphaline,goesline,file=savpath+'combine_goesxray_halpha_flare_lists-res1.sav'

;fix a couple of the columns that have odd blank values

whclassmiss=where(FLRCOMBSTR.HCLASS eq '"' or FLRCOMBSTR.HCLASS eq '""')
FLRCOMBSTR[whclassmiss].HCLASS='Z'

whnoaamiss=where(FLRCOMBSTR.HNOAA eq -1)
FLRCOMBSTR[whnoaamiss].HNOAA=0

write_data_file,FLRCOMBSTR,file=savpath+'combine_goesxray_halpha_flare_lists.csv', $
	header=['#GOES+Halpha+NAR lists are combined to create the most complete flare list with locations possible.','#GOES events are used as a base.','#Missing Values: HST=-1, HPK=-1, HEN=-1, HCLASS=Z, GNOAA=0, HNOAA=0, LAT=1000, LON=1000, HALPHAMATCH=-1, NARMATCH=-1']

;publish the list....

skiptores1:
restore,savpath+'combine_goesxray_halpha_flare_lists-res1.sav',/ver









stop















end