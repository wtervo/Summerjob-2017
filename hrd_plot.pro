;Made by Oskari Tervo as a summer trainee for prof. Jürgen Schmidt
;University of Oulu, Finland, summer 2017
;
;
;data files from https://sbn.psi.edu/archive/cohrd/
;Economou, T. and DiDonna, P., Cassini High Rate Detector V17.0. CO-D-HRD-3-COHRD-V17.0. NASA Planetary Data System, 2017.


function dffrntl, time, events
		
	;function that differentiates particles as a function of time

	;dummy value for the array, calculate number of datapoints
	diff=1d0
	n=n_elements(time)
	
	for i=0,n-1 do begin
		;1st value is at a boundary condition so it's derivative must be handled differently
		if (i eq 0) then begin
			f0=events(i)
			f1=events(i+1)
			f2=events(i+2)
			d0=time(i+1)-time(i+0)
			d1=time(i+2)-time(i+1)

			deri=(((d0+d1)^2)*(f1-f0)-(d0^2)*(f2-f0))/(d0*d1*(d0+d1))
			diff=[diff, deri]
		
		;last value is also handled similarly
		endif else if (i eq n-1) then begin
			f0=events(i)
			f1=events(i-1)
			f2=events(i-2)
			d0=time(i)-time(i-1)
			d1=time(i-1)-time(i-2)

			deri=(((d0+d1)^2)*(f0-f1)-(d0^2)*(f0-f2))/(d0*d1*(d0+d1))
			diff=[diff, deri]
		
		;all the other values go through a "typical" process
		endif else begin
			fi=events(i)
			fip=events(i+1)
			fim=events(i-1)
			dp=time(i+1)-time(i)
			dm=time(i)-time(i-1)

			deri=fi*((dp-dm)/(dp*dm))+(fip*(dm^2)-fim*(dp^2))/(dp*dm*(dp+dm))
			diff=[diff, deri]
		endelse
	
	endfor
	
	;remove the first dummy value		
	diff=diff(1:*)

;return the differentiated list
return, diff
end

function p_corr, detec, et_tme, area

	;function for applying HRD pointing correction to the particle rates detected by HRD
	;1st argument is an array of the rates detected by HRD
	;2nd argument is an et-time array
	;3rd argument is the area of the HRD disc in question

	;structures for data of spacecraft and target moon
	spacecraft={name:'CASSINI',id:-82,gm:-1,reff:-1d0,radii:-1*[1,1,1]}
	target={name:'ENCELADUS',id:-1,gm:-1,reff:-1d0,radii:-1*[1,1,1]}
	;this gives the spice id (integer number) if the body name is provided
	cspice_bodn2c,target.name,id,found
	target.id=id

	;reference frame
	iref='IAU_ENCELADUS'

	;aberration correction is disabled
	abcorr='NONE'

	n=n_elements(et_tme)
	
	;CDA boresight unit vector in cassini_cda-frame
	xyz=[0d0, 0d0, 1d0]

	corr=0d0

	for i=1L, n-1 do begin

		;get v_vector component values from SPICE
		CSPICE_SPKEZ,spacecraft.id,et_tme(i),iref,abcorr,target.id, $
			vec_r_and_v,lt_corr

		;SPICE values in km/s by default - change to m/s
		u=vec_r_and_v(3)*1000d0
		v=vec_r_and_v(4)*1000d0
		w=vec_r_and_v(5)*1000d0
		
		;get transformation vector from cassini_cda frame to Enceladus' frame
		cspice_pxform, iref, 'cassini_cda', et_tme(i), rtt
		
		;transform the boresight unit vector into Enceladus frame
		cspice_mxv, rtt, xyz, bores

		;x, y, and z components of the position vector
		x=bores(0)
		y=bores(1)
		z=bores(2)

		;compute pointing corrected particle value with rate detected by HRD,
		;HRD surface area and dot product of boresight and velocity vectors
		val=detec(i)/(area*(x*u+y*v+z*w))
		corr=[corr, val]
			
	endfor

	;trim array
	corr=corr(1:*)

return, corr
end

function draw_graph, filename, win, day, pl_num, selec

	;function for reading HRD datafiles
	;1st argument is the HDR-datafile
	;2nd argument is the width of the time window from both before and after the closest approach (eq. value of 10 plots events 10mins before and after c/a)
	;3rd argument is the date of the Enceladus flyby in UTC format "YEAR-MONTH-DAY"
	;4th argument is how many of the M1-M4 and m1-m4 particles one wishes to plot, eg. value "2" will plot data for M1, M2, m1, and m2 particles
	;5th argument is "f", "d", or "n" depending on what the user wishes to plot

	device, decomposed=0
	tek_color

	;calls the flybys procedure which defines c/a, v-vector etc. through SPICE
	date=day
	flybys, tra, ca, date

	;as events are only counted every second at most, the 10^-3 seconds are stripped to avoid problems
	;also formats SPICE's UTC day of year format to similar form as in the data files
	ca=strsplit(ca.utc,'.',/extract)
	ca=ca(0)
	ca=strsplit(ca,' // ',/extract)
	ca=ca(0)+"T"+ca(1)

	;opens the file and gets free lun
	file=filename
	openr,lun,/get_lun,file
	line=''
	;dummy values for arrays
	tt=1d0
	b1tot=0d0
	b2tot=0d0
	b3tot=0d0
	b4tot=0d0
	s1tot=0d0
	s2tot=0d0
	s3tot=0d0
	s4tot=0d0
	et_time=0d0
	
	;reads every line of the file in a while-loop
	while not eof(lun) do begin

		readf, lun, line
		temp=strsplit(line,' ',/extract)
		;checks the utc time values of events and defines c/a time with event timer seconds when it finds the the correct value
		utc=temp(5)
		utc=strsplit(utc,'.',/extract)
		utc=utc(0)
		if (utc eq ca) then ca=temp(6)
		;reads the time value of the current line
		tme=temp(6)
		;reads event values for large HRD M1,M2,M3,M4 particles
		b1val=total(fix(temp(9)))
		b2val=total(fix(temp(10)))
		b3val=total(fix(temp(11)))
		b4val=total(fix(temp(12)))
		;reads event values for small HRD m1,m2,m3,m4 particles
		s1val=total(fix(temp(13)))
		s2val=total(fix(temp(14)))
		s3val=total(fix(temp(15)))
		s4val=total(fix(temp(16)))
		;event counter data with the value 256 can cause problems with times especially when differentiating, so these lines are ignored completely
		if (temp(1) eq 256) then continue
		;if current time value is same as the previous one, adds the lhrd/shrd values to the latest value of the respective array instead of creating a new entry
		stop
		if (tme eq tt(-1)) then begin
			
			b1tot(-1)=b1tot(-1)+b1val
			b2tot(-1)=b2tot(-1)+b2val
			b3tot(-1)=b3tot(-1)+b3val
			b4tot(-1)=b4tot(-1)+b4val

			s1tot(-1)=s1tot(-1)+s1val
			s2tot(-1)=s2tot(-1)+s2val
			s3tot(-1)=s3tot(-1)+s3val
			s4tot(-1)=s4tot(-1)+s4val

		;in other cases adds the values incrementally to their appropriate arrays
		endif else begin

			tt=[tt,tme]
			
			b1val=b1tot(-1)+b1val
			b1tot=[b1tot,b1val]
			b2val=b2tot(-1)+b2val
			b2tot=[b2tot,b2val]
			b3val=b3tot(-1)+b3val
			b3tot=[b3tot,b3val]
			b4val=b4tot(-1)+b4val
			b4tot=[b4tot,b4val]

			s1val=s1tot(-1)+s1val
			s1tot=[s1tot,s1val]
			s2val=s2tot(-1)+s2val
			s2tot=[s2tot,s2val]
			s3val=s3tot(-1)+s3val
			s3tot=[s3tot,s3val]
			s4val=s4tot(-1)+s4val
			s4tot=[s4tot,s4val]

			utc=temp(5)
			cspice_str2et, utc, et
			et_time=[et_time, et]

		endelse		
		
	endwhile

	;strips the dummy values and substracts the time value in the beginning from every other time value
	time0=tt(1)
	time=tt(1:*)-time0
	et_time=et_time(1:*)
	
	num=pl_num

	;detector surface areas in m²
	a1=0.0050
	a2=0.0010

	;strips the dummy values of the incremental functions
	b1tot=b1tot(1:*)
	s1tot=s1tot(1:*)

	diff_b1=dffrntl(time, b1tot)
	diff_s1=dffrntl(time, s1tot)

	;if "f" is selected, plots detected particles as a function of time
	if (selec eq "f") then begin

		LHRD_M1=b1tot
		SHRD_M1=s1tot

	;if "d" is selected plots the differentiated time-particle-function
	endif else if (selec eq "d") then begin

		LHRD_M1=diff_b1
		SHRD_m1=diff_s1

	;if  "n" is selected plots number density
	endif else if (selec eq "n") then begin
		
		LHRD_M1=p_corr(diff_b1, et_time, a1)
		SHRD_m1=p_corr(diff_s1, et_time, a2)

	endif

	;same process for M2 and m2 particles
	if (num ge 2) then begin

		b2tot=b2tot(1:*)
		s2tot=s2tot(1:*)

		diff_b2=dffrntl(time, b2tot)
		diff_s2=dffrntl(time, s2tot)

		if (selec eq "f") then begin

			LHRD_M2=b2tot
			SHRD_M2=s2tot

		endif else if (selec eq "d") then begin

			LHRD_M2=diff_b2
			SHRD_m2=diff_s2

		endif else if (selec eq "n") then begin
		
			LHRD_M2=p_corr(diff_b2, et_time, a1)
			SHRD_m2=p_corr(diff_s2, et_time, a2)

		endif

	endif

	;for M3 and m3 particles
	if (num ge 3) then begin

		b3tot=b3tot(1:*)
		s3tot=s3tot(1:*)

		diff_b3=dffrntl(time, b3tot)
		diff_s3=dffrntl(time, s3tot)

		if (selec eq "f") then begin

			LHRD_M3=b3tot
			SHRD_M3=s3tot

		endif else if (selec eq "d") then begin

			LHRD_M3=diff_b3
			SHRD_m3=diff_s3

		endif else if (selec eq "n") then begin
		
			LHRD_M3=p_corr(diff_b3, et_time, a1)
			SHRD_m3=p_corr(diff_s3, et_time, a2)

		endif
	
	endif
	
	;for M4 and m4 particles
	if (num eq 4) then begin

		b4tot=b4tot(1:*)
		s4tot=s4tot(1:*)

		diff_b4=dffrntl(time, b4tot)
		diff_s4=dffrntl(time, s4tot)
		
		if (selec eq "f") then begin

			LHRD_M4=b4tot
			SHRD_M4=s4tot

		endif else if (selec eq "d") then begin

			LHRD_M4=diff_b4
			SHRD_M4=diff_s4

		endif else if (selec eq "n") then begin
		
			LHRD_M4=p_corr(diff_b4, et_time, a1)
			SHRD_m4=p_corr(diff_s4, et_time, a2)

		endif

	endif

	;time of the closest approach
	ca=ca-time0

	;modifying time array to better suit plotting
	tme=(time-ca)/60

	;window layout and plot font size
	!p.multi=[0,1,num]
	!p.charsize=2.5
	
	;time window
	wndw=where(abs(tme) lt win)

	;smoothing for the differentiated function using Savitzky-Golay filter
	smooth=savgol(9, 9, 0, 1)

	;plotting data depending on the number chosen as an argument
	window, 1, xs=1440, ys=900

	for i=1,num do begin

		if (i eq 1) then begin
			ttl=" M1"
			dat=LHRD_M1
		endif else if (i eq 2) then begin
			ttl=" M2"
			dat=LHRD_M2
		endif else if (i eq 3) then begin
			ttl=" M3"
			dat=LHRD_M3
		endif else if (i eq 4) then begin
			ttl=" M4"
			dat=LHRD_M4
		endif
	
		title="Large HRD"+ttl+" particles"
		if (selec eq "d") then title=title+" (differential)
		if (selec eq "n") then title=title+" (number density)

		;plotting of the large HRD data into one window
		bindat, 0.25, tme(wndw), xbin, dat(wndw), ybin

		plot, tme(wndw), convol(dat, smooth), title=title, xtitle="Time from closest approach [minutes]", ytitle="Particles detected"
		oplot, xbin, ybin, psym=10, col=2
		oplot, [0,0], [0,10000000], lin=5

	endfor

	window, 2, xs=1440, ys=900

	for i=1,num do begin

		if (i eq 1) then begin
			ttl=" m1"
			dat=SHRD_m1
		endif else if (i eq 2) then begin
			ttl=" m2"
			dat=SHRD_m2
		endif else if (i eq 3) then begin
			ttl=" m3"
			dat=SHRD_m3
		endif else if (i eq 4) then begin
			ttl=" m4"
			dat=SHRD_m4
		endif

		title="Small HRD"+ttl+" particles"
		if (selec eq "d") then title=title+" (differential)
		if (selec eq "n") then title=title+" (number density)

		;plotting of the small HRD data into a separate window
		bindat, 0.25, tme(wndw), xbin, dat(wndw), ybin

		plot, tme(wndw), convol(dat, smooth), title=title, xtitle="Time from closest approach [minutes]", ytitle="Particles detected"
		oplot, xbin, ybin, psym=10, col=2
		oplot, [0,0], [0,10000000], lin=5

	endfor

	;return window layout to normal
	!p.multi=0

	;close file and free lun
	close,lun
	free_lun,lun	

end

pro hrd_plot

	
	data=draw_graph("hrd_2010_350_365_prc.tab", 10, '2010-12-21', 2, "d")
	;data=draw_graph("hrd_2005_190_199_prc.tab", 10, '2005-07-14', 2, "d")
	

end
