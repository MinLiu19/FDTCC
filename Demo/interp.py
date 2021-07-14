#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    A simple scirpt to interp data in frequency domain

Created on Sun May 16 16:00:13 2021
@author: mliu
"""
import numpy as np
from obspy import read, UTCDateTime
from scipy.fftpack import fft,ifft
import os
import sys

catalog = "./hypoDD.pha"
data = "./waveforms/"
delta1 = float(sys.argv[1])
delta2 = float(sys.argv[2])
time = 50
freqmin = 1
freqmax = 15
npts1 = time/delta1
npts2 = time/delta2
with open(catalog,"r") as ots:
    for ot in ots:
        if(ot.split()[0]=="#"):
           jk1, year, month, day, hour, minute, tmp, evla, evlo, evdp, jk2, jk3, jk4, jk5, eventid  = ot.split()
           print(year+"/"+month+"/"+day+" "+hour+":"+minute+":"+tmp)
           sec,msec=tmp.split(".")
           outdir = data + str(eventid)
           indir = data + str(year) + str(month) + str(day)
           if not os.path.exists(outdir):
               os.mkdir(outdir)
           starttime = UTCDateTime(int(year),int(month),int(day),int(hour),int(minute),int(sec),int(msec+"000"))
           endtime = starttime + time
           for staname in os.listdir(indir):
               try:
                   st = read(indir + "/" + staname,starttime=starttime,endtime=endtime)
                   st.detrend("demean")
                   st.detrend("linear")
                   st.filter('bandpass',freqmin=freqmin,freqmax=freqmax,zerophase=True)
                   st_fft = fft(st[0].data[0:int(npts1)])
                   zero = np.full(int(npts2-npts1),np.complex64(0))
                   st_zero = st_fft[0:int(npts1/2)]
                   st_zero = np.hstack((st_zero,zero))
                   st_zero = np.hstack((st_zero,st_fft[int(npts1/2):int(npts1)]))
                   st[0].data = ifft(st_zero)
                   st[0].stats.delta = delta2
                   st[0].write(outdir + "/"+st[0].stats.network+"."+st[0].stats.station+"."+st[0].stats.channel,format="SAC")
               except:
                   pass
