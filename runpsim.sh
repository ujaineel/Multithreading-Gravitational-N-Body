#!/bin/bash
#Workflow script for particle sim
mydir=$PWD
rm anim.gif
rm -f timedat.*

#execute code
./nBody $1 $2

#performn visualization
test=0;
tops =15;
while (($test <= $2-1))
do
echo "set xrange [0:1]
set title \"$1 particles at timestep $test\"
set yrange [0:1]
set timestamp top $tops
set grid
set term gif size 540,700
set output '$test.gif'
plot \"timedat.0\" i $test u 4:5 pt 3 ps 1 t \"Node\;" >data_$test.gnu
gnuplot data_$test.gnu
let test=$test+1
done


#cleanup
rm -f *.gnu
ls *.gif | sort -nk1 | xargs ./gifmerge -10 -l0 >anim.vid
rm -f *.gif
mv anim.vid anim.gif
rm -f timedat.*

