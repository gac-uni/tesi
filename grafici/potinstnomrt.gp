set terminal epslatex font "Times-New-Roman, 12"
# set terminal qt persist
set out "potinstnomrt.tex"

unset border
unset key
unset grid

set xtics axis
set ytics axis

set style arrow 1 head filled size char 1,20,50 lw 1.5

set arrow 1 from 0,0 to 27000,0 arrowstyle 1
set arrow 2 from 0,-100 to 0, 450 arrowstyle 1

# set xrange[-13:13]
set yrange[-100:400]

set tmargin 4
set rmargin 8
set bmargin 5

set ytics scale 0.5

set xtics ("$1\\cdot10^5$" 5000,"$2\\cdot10^5$" 10000,"$3\\cdot10^5$" 15000 ,"$4\\cdot10^5$" 20000 ,"$5\\cdot10^5$" 25000 ) 
set xtics offset 0, -1

set ylabel '$(1/N_{p})V_{lj}(\vb*{r}_i) [K]$'
set xlabel 'Configurations in the sample'
set xlabel offset 0,-1

plot "< awk 'BEGIN{i=0}{if(i%20==0){print $0}; i++}' ../varmcnomrt/mt_val_dim/pot_10.dat" u 0:1 w p pt 7 ps 0.1 lc rgb "#0000FF" 
 