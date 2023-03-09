set terminal epslatex font "Times-New-Roman, 12"
# set terminal qt persist
set out "potljpresentazione.tex"
set samples 1000

unset border
set key
set key at screen 1, 0.9
unset grid

set xtics axis
set ytics axis

set style arrow 1 head filled size char 1,20,50

set arrow 1 from graph 0,0.5 to screen 0.99, graph 0.5 arrowstyle 1
set arrow 2 from graph 0,0 to graph 0, screen 0.99 arrowstyle 1

eps=10.22
sigma=2.556
f(x) = 4*eps*((sigma/x)**12-(sigma/x)**6)

set tmargin 3
set rmargin 8

set xrange [0:8]
set yrange [-20:20]

set xlabel "$r[\\text{\\AA}]$"
set ylabel "$V_{lj}[\\textrm{K}]$" #offset 8, 0, 0

plot f(x) lc rgb "#0000FF" t "$V_{lj} = 4\\varepsilon \\qty[ \\qty(\\frac{\\sigma}{r})^{12}-\\qty(\\frac{\\sigma}{r})^{6}]$"