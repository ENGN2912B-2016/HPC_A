set pm3d 
set palette defined (10 'black', 30 'gray')
set terminal gif animate
set output "ternary_thread_equal_ie.gif"
unset border
unset xtics
unset ytics
unset ztics
unset colorbox
n=100000
i=0
load "animate3.gp"
set output
