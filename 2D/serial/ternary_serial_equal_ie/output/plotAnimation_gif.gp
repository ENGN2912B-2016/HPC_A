set cbrange [10:30]
set size square
set palette defined (10 'black', 30 'gray')
set terminal gif animate
set output "ternary_serial_equal_ie.gif"
unset border
unset xtics
unset ytics
unset colorbox
n=100000
i=0
load "animate2.gp"
set output

