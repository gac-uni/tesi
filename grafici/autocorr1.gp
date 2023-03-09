set terminal epslatex font "Times-New-Roman, 12"

set out "autocorr1.tex"

f(x) = exp(-x)

set xrange [0:10]
unset border
unset key 
unset grid

set tmargin 3
set rmargin 8

set xtics axis
set ytics axis
set tics in 
set tics scale 1
set mxtics 20
set mytics 10
set ytics 1 

set style arrow 1 head filled size char 1,20,50
set arrow 1 from graph 0,0 to screen 0.99, graph 0 arrowstyle 1
set arrow 2 from graph 0,0 to graph 0, screen 0.99 arrowstyle 1

set label "\\huge $\\tilde{C}_k = e^{-\\frac{k}{\\tau}}$" at graph 0.6, 0.7 font "Times,80"
set xlabel "$k$"
set ylabel "$\\tilde{C}_k$"



plot f(x) lc rgb "#0000FF"

