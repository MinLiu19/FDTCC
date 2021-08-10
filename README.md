			Fast double-difference cross-correlation (FDTCC)
							Authors: Min Liu & Miao Zhang
							m.liu@dal.ca & miao.zhang@dal.ca
	1.Usage
		FDTCC -C(ife/ifd/ifp) -W(wb/wa/wf/wbs/was/wfs) -D(delta/threshold/thre_SNR/thre_shift) -G(trx/trh/tdx/tdh) 
		-B(low/high) -F(f)
   	--------------------------------------explanation----------------------------------------
		-C: specify the path of event.sel, dt.ct and phase.dat (1: yes, 0: default names)
		-W: waveform window length before and after picks and their maximum shift length
		-D: sampling interval, CC threshold, SNR threshold, maximum abs(t1-t2) of the two picks
		-G: ranges and grids in horizontal direction and depth (in traveltime table)
		-F: input data format (-3: continuous data; -5: event segments)
		-B: waveform bandpass filtering (e.g., 2/8; -1/-1: no filter applied).
      	   	SAC name format: date/net.sta.comp, e.g., 20210101/AA.BBBB.HHZ
                           or eventID/net.sta.comop, e.g., 8/AA.BBBB.HHZ).
		staDir: stations directory
       		tttDir: travel-time directory
       		wavDir: waveform directory
       		eveDir: event.sel (optional)
       		dctDir: dt.ct (optional)       
		paDir:	phast.dat (optional)

	2.Input file
		phase.dat		# in hypoDD format
		event.sel 		# in hypoDD format
		dt.ct     		# in hypoDD format
		station.dat		# in REAL format
		waveform data		# raw data
		travel-time table	# in REAL format

	3.Outfile
		dt.cc			# can be used in hypoDD, Growclust and tomoDD directly

	4.Demo
		$bash runFCC.sh
![image](https://github.com/MinLiu19/FDTCC/blob/main/Workflow.jpg)
