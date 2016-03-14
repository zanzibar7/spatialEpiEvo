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

echo $y $x | gnuplot -persist > /dev/null  &
#echo $y $x

exit

#gnuplot#
set term x11;
set nokey;
set logscale xy;
plot 0 
#gnuplot#

