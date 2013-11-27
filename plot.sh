#!/bin/sh

#filename="data2"
filename="$1"

echo $filename

y=`sed -n '/^#gnuplot#/,/#gnuplot#/p' $0| grep -v "#" `

let a=0
e=`cat $filename | grep '#' | tail -n 1 | cut -d ' ' -f 3`
let e=e+1
echo $e
while [ $a -lt $e ]; do
	x+=", \"$filename\" using 5:6 index $a w l"
	let a=1+a
done

echo $y $x | gnuplot -persist &
#echo $y $x
exit

#set logscale xy; set yrange [.1:1]; set xrange [.1:10];
#gnuplot#
set term x11;
set nokey;
set yrange [0:1]; set xrange [0:10];
set xlabel "Transmission rate";
set ylabel "Infection period";
set samples 1e3;
plot 1/(8*x) w l, 2/(8*x) w l, 4/(8*x) w l
#gnuplot#

