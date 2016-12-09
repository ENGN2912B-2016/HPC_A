#!/bin/bash

k=0;
echo "set cbrange [10:30]" >> plotAnimation.gp
echo "set size square" >> plotAnimation.gp
echo "set palette defined (10 'black', 30 'gray')" >> plotAnimation.gp
echo "set terminal png" >> plotAnimation.gp
echo "unset border" >> plotAnimation.gp
echo "unset xtics" >> plotAnimation.gp
echo "unset ytics" >> plotAnimation.gp
echo "unset colorbox" >> plotAnimation.gp
while [ $k -lt 100000 ]; do
echo "set output 'timestep$k.png'" >> plotAnimation.gp
echo "plot \"timestep$k.dat\" matrix with image" >> plotAnimation.gp
echo "pause 0.05" >> plotAnimation.gp
let k=$k+1000;
done