set terminal epslatex font "Times-New-Roman, 12"
#set terminal qt persist
set out "blocking.tex"

f(x) = 2/pi*atan(x)

set xrange [0:70]
unset border
unset key 
unset grid

set tmargin 3
set rmargin 8

unset xtics
set xtics format ""
set xtics scale 0

set ytics format ""
set ytics axis
set ytics scale 1
set ytics add ('$\sigma_{true}^2(S_O(X))$' 1) 

set style arrow 1 head filled size char 1,20,50
set style arrow 2 nohead dt "."

set arrow 1 from graph 0,0 to screen 0.99, graph 0 arrowstyle 1
set arrow 2 from graph 0,0 to graph 0, screen 0.99 arrowstyle 1
set arrow 3 from 0,1 to graph 1, 1 arrowstyle 2 

set xlabel "Length of the Blocks $B$"
set ylabel '$\sigma^2(S_O(B))$' offset 8, 0, 0

plot f(x) lc rgb "#0000FF" dt "-"  