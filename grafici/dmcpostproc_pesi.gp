set terminal epslatex font "Times-New-Roman, 12"
# set terminal qt persist
set out "dmcpostproc_pesi.tex"

unset border
unset key
unset grid

unset xtics
set xtics scale 0
set ytics axis

set style arrow 1 head filled size char 1,20,50 lw 1.5
set style line 1 dt 4 lw 2 lc "red"
set style line 2 dt 4 lw 2 lc "blue"
set style line 3 lc rgb "#FEFF0000" 
set style line 4 lc rgb "#FE0000FF"
#set style points 1 pt 7 ps 0.1 lc rgb "#EEFF0000"
#set style points 2 pt 7 ps 0.1 lc rgb "#EE0000FF"

set arrow 1 from 0,1*10**(-60) to 1100,1*10**(-60) arrowstyle 1
set arrow 2 from 0,1*10**(-60) to 0, 10000 arrowstyle 1

# set object 1 circle at 470869,-11 size char .2 fillcolor rgb "red" fillstyle solid
# set object 2 circle at 470869,-13 size char .2 fillcolor rgb "blue" fillstyle solid

set xrange[0:1100]
# set yrange[-15:30]

set label '|' at 200, 1*10**(-60) font "Times, 8"
set label '|' at 400, 1*10**(-60) font "Times, 8"
set label '|' at 600, 1*10**(-60) font "Times, 8"
set label '|' at 800, 1*10**(-60) font "Times, 8"
set label '|' at 1000, 1*10**(-60) font "Times, 8"

# set tmargin 4
# set rmargin 8
# set bmargin 5

set ylabel 'Weights'
set xlabel 'Steps'

set xtics ("$1\\cdot10^5$" 200,"$2\\cdot10^5$" 400,"$3\\cdot10^5$" 600 ,"$4\\cdot10^5$" 800 ,"$5\\cdot10^5$" 1000 )
set logscale y
plot for [i=1:40] "../dmc_postp/outpesi_smol.dat" u 0:i w l