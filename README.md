			Fast double-difference cross-correlation (FDTCC)
							Authors: Min Liu & Miao Zhang
							mliu@cuhk.edu.hk & miao.zhang@dal.ca
	1.Usage
		FDTCC -C(ife/ifd/ifp) -W(wb/wa/wf/wbs/was/wfs) -D(delta/threshold/thre_SNR/thre_shift) -G(trx/trh/tdx/tdh) 
		-B(low/high) -F(f)
   	--------------------------------------explanation----------------------------------------
		-C: specify the path of event.sel, dt.ct and phase.dat (1: yes, 0: default names)
		-W: waveform window length before and after picks and their maximum shift length
		-D: sampling interval, CC threshold, SNR threshold, maximum abs(t1-t2) of the two picks
		-G: ranges and grids in horizontal direction and depth (in traveltime table)
		-F: input data format (0: continuous data; 1: event segments)
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
		
	If you use FDTCC in research that you submit, please cite the two papers:
		1) Multistage Nucleation of the 2021 Yangbi MS 6.4 Earthquake and Its Foreshocks, Journal of Geophysical Research: Solid Earth, 2022, https://doi.org/10.1029/2022JB024091.
		2) Investigation of the 2013 Eryuan, Yunnan, China MS 5.5 Earthquake Sequence: Aftershock Migration, Seismogenic Structure and Hazard Implication, Tectonophysics, 2022, https://doi.org/10.1016/j.tecto.2022.229445
	
	FDTCC is integrated into the earthquake detection and location workflow - LOC-FLOW, more details please see our SRL paper:  
		LOC-FLOW: An End-to-End Machine-Learning-Based High-Precision Earthquake Location Workflow, Seismological Research Letters, 2022, https://doi.org/10.1785/0220220019
![image](https://github.com/MinLiu19/FDTCC/blob/main/Workflow.jpg)
