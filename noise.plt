#!/bin/bash

gnuplot -persist <<EOF

set pm3d map
set size square
set palette grey

set terminal x11 enhanced "Times-Roman, 40"
unset key

set xlabel 'x' font "Times-Roman-Italic, 40"
set ylabel 'y' font "Times-Roman-Italic, 40"

set xtics 250
set ytics 250
set cbtics 5
unset colorbox

splot '$1.dat' u 1:2:3

set terminal post "Times-Roman" 40
set out '$1.ps'
splot '$1.dat' u 1:2:3

EOF