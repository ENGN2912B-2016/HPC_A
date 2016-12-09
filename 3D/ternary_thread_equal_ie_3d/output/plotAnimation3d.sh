#!/bin/bash

k=0;
echo "set pm3d" >> plotAnimation3d.gp
echo "" >> plotAnimation3d.gp
echo "set terminal png" >> plotAnimation3d.gp
echo "unset border" >> plotAnimation3d.gp
echo "unset xtics" >> plotAnimation3d.gp
echo "unset ytics" >> plotAnimation3d.gp
echo "unset ztics" >> plotAnimation3d.gp
echo "unset colorbox" >> plotAnimation3d.gp
while [ $k -lt 100000 ]; do
echo "set output 'timestep$k.png'" >> plotAnimation3d.gp
echo "splot \"timestep$k.dat\" u 1:2:3:4 with points palette notitle" >> plotAnimation3d.gp
echo "pause 0.05" >> plotAnimation3d.gp
let k=$k+1000;
done