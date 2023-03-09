set terminal epslatex font "Times-New-Roman, 12"
# set terminal qt persist
set out "dmcpostproc.tex"

unset border
# unset key
unset grid
set key at 500000, -6

set xtics axis
set ytics axis

set style arrow 1 head filled size char 1,20,50 lw 1.5

set arrow 1 from 0,0 to 520000,0 arrowstyle 1
set arrow 2 from 0,-15 to 0, 18 arrowstyle 1

# set xrange[-13:13]
# set yrange[0:40]

set tmargin 4
set rmargin 8
set bmargin 5

set xtics ("$1\\cdot10^5$" 100000,"$2\\cdot10^5$" 200000,"$3\\cdot10^5$" 300000 ,"$4\\cdot10^5$" 400000 ,"$5\\cdot10^5$" 500000 )

set ylabel 'Energy $[\unit{\kelvin}]$'
set xlabel 'Configurations in the sample'

# set ytics scale 0.5

plot "../dmc_postp/outvaloricum.dat" u 0:1 w l lc "red" t 'Cumulative average of $E_T$', "../dmc_postp/outvaloricum.dat" u 0:2 w l lc "blue" t 'Weighted average of the energy'