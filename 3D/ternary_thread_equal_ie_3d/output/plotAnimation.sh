#!/bin/bash

k=0;
echo "set pm3d" >> plotAnimation3d.gp
echo "set palette defined (10 'black', 30 'gray')" >> plotAnimation3d.gp
echo "set terminal png" >> plotAnimation3d.gp
echo "unset border" >> plotAnimation3d.gp
echo "unset xtics" >> plotAnimation3d.gp
echo "unset ytics" >> plotAnimation3d.gp
echo "unset colorbox" >> plotAnimation3d.gp
while [ $k -lt 100000 ]; do
echo "set output 'timestep$k.png'" >> plotAnimation.gp
echo "plot \"timestep$k.dat\" matrix with image" >> plotAnimation.gp
echo "pause 0.05" >> plotAnimation.gp
let k=$k+1000;
done