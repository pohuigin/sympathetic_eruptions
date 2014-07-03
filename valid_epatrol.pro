;+
;Name:
;
;     valid_epatrol
;     
;Purpose:
;
;     Track events at in a given time range
;
;Calling Sequence:
;
;     valid_epatrol, [ minvalue, maxvalue ]
;
;Optional Inputs:
;     
;     minvalue = The cluster to begin at
;     maxvalue = The cluster to stop at
;
;Outputs:
;
;     Will output values into a csv under the default title of: fe_times_sources.txt
;     
;Some notes:
;     The program will read in values from compile_eventlist_visual_search_clusters_test.txt by default
;     The default cadence is 1 frame every 12 minutes
;     The first cluster is cluster #0
;
;Author:
;     
;     George Lee
;     Intern At LMASAL
;     Email: glee@lmsal.com
;
;
;Running this program:
;     Type in valid_epatrol
;To pick up where you left off,
;     Type in valid_epatrol, cluster
;
;This program will also work when working on clusters in a specified range
;     valid_epatrol, min, max
;-

pro run_epatrol, whichevent, infile, logfilein, verb=verb


;Default Output  
thisline={available:'',cluster:'', flare:'', region:'',tinit:'',xsource:0.0D,ysource:0.0D,xloc1:0.0D,yloc1:0.0D,xloc2:0.0D,yloc2:0.0D}

  ;Sets up read/write variables
  logfile=logfilein
  file = infile
 
  ;Set up identification variables

  clust = file.CLUSTID
  clusterid = clust[whichevent]
  thisline.cluster = clusterid
  
  ;set parameters for movie
  
  cadence=6 ;x2 = 12 minutes between images
  stmovie=0. ; = 0 hours before Eruption Patrol event start time

;Read in eruption patrol events using HEK API

  ;Set up time variables
  
  d1 = file.DATEST
  d2 = file.DATEEN
  g1 = d1[whichevent]
  g2 = d2[whichevent]
  t1=anytim(g1)-stmovie*3600.
  t2=anytim(g2)
  
  ;Directory where program will fetch the images

  dir = '/archive/sdo/AIA/synoptic'

  ;Finds http:// directories for the given files

  ff=getjsoc_synoptic_read(t1,t2,wave = 304,/nodat, $ ; cadence=20*60., $
    info=info,file=files,remfile=remfiles,outind=inds)
   
  ;; Checks to see if there is data for the given date
  
  IF (ff[0] ne '') THEN BEGIN

  ;Sets up variables about number of images to be read in
  
  nimg=n_elements(ff)
  imgnums=findgen(nimg)
  ff=ff[where(imgnums mod cadence eq 0)]
  nimg=n_elements(ff)
  
  ;example return link :   http://jsoc.stanford.edu/data/aia/synoptic//2011/01/02/H1700/AIA20110102_1716_0304.fits
  ;;goal: extract /2011/01/02/H1700/AIA20110102_1716_0304.fits
  ;done by doing: test = 'http://jsoc.stanford.edu/data/aia/synoptic//2011/01/02/H1700/AIA20110102_1716_0304.fits'
  ;c = strmid(teststr, 43)
  
  ;;Converts http:// directories to a /yyyy/mm/dd/Hxxxx/ format
  
  ffloc = ''
  FOR i=0, nimg-1 DO BEGIN
    IF (i eq 0) THEN BEGIN
      ffloc = [dir+ strmid(ff[i], 43)]
    ENDIF ELSE BEGIN
      ffloc = [ffloc, dir+ strmid(ff[i], 43)]
    ENDELSE

  ENDFOR 
  
  ;read in the images

  read_sdo,ffloc,indarr,datarr

  ;Sets color table to one that is more natural looking for the sun

  loadct,3

  ;display initial image, and have user pick out the source


  keyin = ''
  counter = 0

  print, 'Ready!'
  print, 'q - quit, s - skip current cluster, d - back, f - forward, m - play movie, g - input time'
  WHILE (keyin ne 'q') DO BEGIN
    
    ;Sets up array for scrolling through the images
    
    index2map,indarr[counter],datarr[*,*,counter],imap
    plot_map,imap,/log
    time = indarr[counter].date_obs
    thisline.tinit = time
    keyin = get_kbrd()
     if keyword_set(verb) THEN BEGIN
      print, keyin
     endif
     
     ;Plays images in a semi-fast sequence to give the appearance of a movie
     
     IF (keyin eq 'm') THEN BEGIN
;        for i=0, nimg-1,1 do begin & plot_image,alog10(datarr[*,*,i]),title = indarr[i].date_obs & &wait,0.1 & endfor
        for i=0, nimg-1,1 do begin
          index2map,indarr[i],datarr[*,*,i],imap
          plot_map,imap,/log
          wait,0.1
         endfor

     ENDIF
     
    ;Takes in Event Related Statistics

    IF (keyin eq 'g') THEN BEGIN
      print,'Click on the source location of the eruption, Right click to confirm your selection'
      xstore = 0
      ystore = 0
      mousenum = 0


;    ;Picks Location of Event
      print, 'Plot 1'
      WHILE (mousenum ne 4) DO BEGIN
        cursor,xinit,yinit
        mousenum = !mouse.button

        IF(mousenum eq 1) THEN BEGIN
          xstore = xinit & ystore = yinit
          ERASE
           plot_map,imap,/log
          vline,xinit & hline,yinit
        ENDIF
        IF(mousenum eq 4) THEN BEGIN
          xinit = xstore & yinit = ystore
          mousenum = 0
          mousenum2 = 0
          BREAK
        ENDIF
        wait,0.05
      ENDWHILE
      vline,xinit & hline,yinit
      print,xinit,yinit
      wait, 1.0
      thisline.xsource = xinit
      thisline.ysource = yinit
    

   ;Creates bounding box of the event
   ERASE
   plot_map,imap,/log
       print, 'Click two times to select a bounding box for the event, Right click to confirm your selection'
       print, 'Plot 2'
       mousenum2 = 0
       xstore1 = 0
       ystore1 = 0
       xstore2 = 0
       ystore2 = 0
       counter1 = 0
       counter2 = 0
       
;Finds bounding box 

      WHILE(mousenum2 ne 4) DO BEGIN
          ;crosshair 1
        
          ;find point for 1
          cursor, xstore1, ystore1
          mousenum2 = !mouse.button
          IF (mousenum2 eq 1) THEN BEGIN
          IF(counter1 eq 0) THEN BEGIN
            vline,xstore1 & hline,ystore1
            print, '1 done'
            print, xstore1 & print, ystore1
            wait,0.5        
            counter1 += 1
          ENDIF ELSE BEGIN
          ;Replots points
            ERASE
            plot_map,imap,/log
            vline,xstore1 & hline,ystore1
            vline,xstore2 & hline,ystore2
            print, '1 done'
            print, xstore1 & print, ystore1
            wait,0.5
          ENDELSE
          ENDIF
      
          ;crosshair 2
         
          ;find point for 2
         cursor, xstore2, ystore2
         mousenum2 = !mouse.button
          ;crosshair 2
          IF (mousenum2 eq 1) THEN BEGIN
          IF(counter2 eq 0) THEN BEGIN
            vline,xstore2 & hline,ystore2
            print, '2 done'
            print, xstore2 & print, ystore2
            wait,0.5
            counter2 +=1
          ENDIF ELSE BEGIN
          ;Replots points
            ERASE
            plot_map,imap,/log
            vline,xstore1 & hline,ystore1
            vline,xstore2 & hline,ystore2
            print, '2 done'
            print, xstore2 & print, ystore2
            wait,0.5
          ENDELSE
          ENDIF
          
      ENDWHILE
      
      
      thisline.xloc1 = xstore1 & thisline.yloc1 = ystore1
      thisline.xloc2 = xstore2 & thisline.yloc2 = ystore2
           
;      store1 = [xstore1,ystore1]
;      store2 = [xstore2,ystore2]
;      loc = [store1,store2]
      
      ERASE
      plot_map,imap,/log
      
      while (0 ne 1) DO BEGIN
        print, 'Flare or no flare? (f/n)'
      type = get_kbrd()
      IF (type eq 'f') or (type eq 'n')THEN BEGIN
        thisline.flare = type
        BREAK
      ENDIF ELSE BEGIN
        print, 'Input Is Not Valid'
      ENDELSE
      endwhile
      
      while (0 ne 1) DO BEGIN
        print, 'Event from an active region or quiet sun? (a/q)'
        type = get_kbrd()
        IF (type eq 'a') or (type eq 'q')THEN BEGIN
          thisline.region = type
          BREAK
        ENDIF ELSE BEGIN
          print, 'Input Is Not Valid'
        ENDELSE
      endwhile
      
;Writes into data file
      thisline.available = 'y'
      write_data_file,thisline,filecsv=logfile,/append
      print, 'Written!'
    ENDIF



    IF (keyin eq 'f') THEN BEGIN
      IF (counter NE nimg-1) THEN counter +=1
    ENDIF

    IF (keyin eq 'd') THEN BEGIN
      IF (counter GT 0) THEN counter -=1
    ENDIF

    ;Skips current cluster

    IF (keyin eq 's') THEN BEGIN
      thisline.available = 's'
      thisline.tinit = g1
      write_data_file,thisline,filecsv=logfile,/append
      print, 'Written!'
      BREAK
    ENDIF

  ENDWHILE

  
  
  ;;Prints out log file for dates with no data
  ENDIF ELSE BEGIN
    print, 'There is no data for the date specified.'
     thisline.available = 'n'
    thisline.tinit = g1
    write_data_file,thisline,filecsv=logfile,/append
    print, 'Written!'
  ENDELSE


  print, 'Finished!'
  stop

end

pro valid_epatrol, minvalue, maxvalue

  logfile='fe_times_sources.txt'
  infile= 'compile_eventlist_visual_search_clusters_test.txt'

  file = read_data_file(infile)
  
  ;Will automatically create a blank logfile in the current IDL directory if there is no existing one
  if(~file_test(logfile)) then begin
    thisline={available:'',cluster:'', flare:'', region:'',tinit:'', xsource:0.0D, ysource:0.0D, xloc1:0.0D, yloc1:0.0D, xloc2:0.0D, yloc2:0.0D}
    write_data_file,thisline,filecsv=logfile,/nodata
    CD, Current=theDirectory
    print, 'New File Written To ' + theDirectory
  endif

;Automatically sets minimum value to 0 if not specified
  if (minvalue eq !NULL) then begin
    minval = 0
  endif else begin
    minval = minvalue
  endelse
;Automatically sets the maximum value to the last value in the file read in
  if (maxvalue eq !NULL) then begin
    maxval = n_elements(file)
  endif else begin
    maxval = maxvalue
  endelse

;
  for i= minval,maxval,1 do begin
    run_epatrol, i, file,logfile
    clustnum = STRTRIM(i, 2)
    print,'Cluster #'+ clustnum +' is complete'
  endfor
end