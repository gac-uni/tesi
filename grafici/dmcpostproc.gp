set terminal epslatex font "Times-New-Roman, 12"
# set terminal pdf font "Times-New-Roman, 12"
# set terminal pngcairo font "Times-New-Roman, 12"
# set terminal qt persist
set out "dmcpostprocpresentazione.tex"

unset border
# unset key
unset grid
set key at 850000, 30

set xtics axis
set ytics axis

# set xtics scale 0
# unset xtics
# set ytics scale 0
# unset ytics

set style arrow 1 head filled size char 1,20,50 lw 1.5
set style line 1 lw 1.5 lc "green"
set style line 2 lw 1.5 dt 3 lc "blue"
# set style point 1 pt 7 ps 0.1 lc rgb "#5BFF93"
# set style point 2 pt 7 ps 0.1 lc rgb "##02D8DF"

set arrow 1 from 0,0 to 520000,0 arrowstyle 1
set arrow 2 from 0,-15 to 0, 35 arrowstyle 1

# set object 1 circle at 476500,-11 size char .2 fillcolor rgb "red" fillstyle solid
# set object 2 circle at 476500,-13 size char .2 fillcolor rgb "blue" fillstyle solid

# set xrange[-13:13]
set yrange[-15:30]

set tmargin 4
set rmargin 20
set bmargin 5
set tics front

set xtics ('$1 \cdot 10^5$' 100000,'$2 \cdot 10^5$' 200000,'$3 \cdot 10^5$' 300000 ,'$4 \cdot 10^5$' 400000 ,'$5 \cdot 10^5$' 500000 )

set ylabel 'Energy $[ \unit{\kelvin} ]$'
set xlabel 'Configurations in the sample'

# set ytics scale 0.5

plot "< awk 'BEGIN{i=0}{if(i%50==0){print $0}; i++}' ../dmc_postp/outvalori.dat" u ($0*50):1 w p pt 7 ps 0.1 lc rgb "#5BFF93" t '$E_T$',\
    "< awk 'BEGIN{i=0}{if(i%50==0){print $0}; i++}' ../dmc_postp/outvalori.dat" u ($0*50):2 w p pt 7 ps 0.1 lc rgb "#02D8DF" t '$\ev{K(r_i)}+\ev{V_{lj}(r_i)}$',\
    "../dmc_postp/outvaloricum.dat" u 0:1 w l ls 1 t 'Cum. avg of $E_T$',\
    "../dmc_postp/outvaloricum.dat" u 0:1 w l ls 2 t 'Cum. avg of $\ev{K(r_i)}+\ev{V_{lj}(r_i)}$'
