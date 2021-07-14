#!/bin/bash
# create the dt.cc file from continuous SAC files or event segments

#-----------------------------Parameters setting-----------------------------
#waveform window length before and after trigger
#
#  -------|-----------|------------------|------
#             (wb)   pick      (wa)
# For P phase CC, if your "wa" is larger than 0.9*(ts-tp), it will be replaced
# by 0.9*(ts-tp) to make sure you don't include S phase.
# For S phase CC, if your "wb" is larger than 0.5*(ts-tp), it will be replaced
# by 0.5*(ts-tp) to make sure you don't include P phase.
#waveform window length before and after picks and their maximum shift length
W=0.2/1.0/0.3/0.5/1.5/0.5
#sampling interval, CC threshold, SNR threshold, maximum abs(t1-t2) of the two picks
D=0.01/0.7/1/2
#ranges and grids in horizontal direction and depth (in traveltime table)
G=3/20/0.02/2
#specify the path of event.sel, dt.ct and phase.dat (1: yes, 0: default names)
C=1/1/1
#input data format (-3: continuous data -5: interpolated event segments)
F=-5
if [ $F == -5 ];then
echo "cut and interpolate event segments"
delta=`echo $D|gawk -F/ '{print $1}'`
python interp.py 0.01 $delta #interp data from 0.01 to $delta
fi
#BP filter, low and high B=-1/-1 will not fiter the data 
B=2/8

staDir=./station.dat
tttDir=./tt_db/ttdb.txt
wavDir=./waveforms
eveDir=./event.sel
dctDir=./dt.ct
paDir=./hypoDD.pha
../bin/FDTCC -F$F -B$B -C$C -W$W -D$D -G$G $staDir $tttDir $wavDir $eveDir $dctDir $paDir
echo ../bin/FDTCC -F$F -B$B -C$C -W$W -D$D -G$G $staDir $tttDir $wavDir $eveDir $dctDir $paDir
rm Input*
