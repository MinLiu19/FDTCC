#!/bin/bash
Trace=/Users/mliu/2021/software/FDTCC/Demo/waveforms/20190704
Template=/Users/mliu/2021/software/FDTCC/Demo/waveforms
for event in `cat ./catalog.dat|gawk '{if($1==1||$1==7)print substr($2,1,4)"/"substr($2,6,2)"/"substr($2,9,2)"/"substr($3,1,2)"/"substr($3,4,2)"/"substr($3,7)"/"$4"/"$5"/"$6"/"$7"/"$1}'`
do
	echo $event
	year=`echo $event | gawk -F/ '{print $1}'`
     	month=`echo $event | gawk -F/ '{print $2}'`
     	day=`echo $event | gawk -F/ '{print $3}'`
     	hour=`echo $event |gawk -F/ '{print $4}'`
     	min=`echo $event | gawk -F/ '{print $5}'`
     	tmp=`echo $event |gawk -F/ '{print $6}'`
     	sec=`echo $tmp |gawk -F. '{print $1}'`
	msec=`echo $tmp |gawk -F. '{print $2}'`	
	elat=`echo $event |gawk -F/ '{print $7}'`
	elon=`echo $event |gawk -F/ '{print $8}'`
	edep=`echo $event |gawk -F/ '{print $9}'`
	emag=`echo $event |gawk -F/ '{print $10}'`
	ID=`echo $event |gawk -F/ '{print "test_"$11}'`
cd $Trace
#mkdir -p $Template/$year$month$day$hour$min$sec
mkdir -p $Template/$ID
#ID=$year$month$day$hour$min$sec
#mkdir -p $Template/$ID
mins=`echo $min|gawk '{print $1}'`
for sta in `ls` 
do
sac<<EOF
r $sta
ch o gmt $year 185 $hour $mins $sec $msec
setbb temp ( MULTIPLY &1,o -1 )
ch allt %temp
w $Template/$ID/$sta.T1
cut 0 50
r $Template/$ID/$sta.T1
ch b 0
ch evla $elat evlo $elon evdp $edep user0 $emag 
rmean;rtr;taper
bp c 1 15 p 2 n 4 
w $Template/$ID/${sta}
quit
EOF
rm $Template/$ID/*T1
done
done
