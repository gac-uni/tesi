set terminal epslatex font "Times-New-Roman, 12"
# set terminal qt persist

set out "birthdeatden.tex"

unset border
# unset key
unset grid
set key 

set xtics axis
set ytics axis

set tmargin 5
set rmargin 5

set style arrow 1 head filled size char 1,20,50 lw 1.5

set arrow 1 from 0,0 to 520000,0 arrowstyle 1
set arrow 2 from 0,-10 to 0, 12 arrowstyle 1

set yrange [-10:10]

set xtics ("$1\\cdot10^5$" 100000,"$2\\cdot10^5$" 200000,"$3\\cdot10^5$" 300000 ,"$4\\cdot10^5$" 400000 ,"$5\\cdot10^5$" 500000 ) 

set ylabel 'Energy [ K ]'
set xlabel 'Configurations in the sample'

plot "< awk 'BEGIN{i=0}{if(i%20==0){print $0};i++}' ../birthdeath/valen.dat" u ($0*20):($1*10.22/64) w p pt 7 ps 0.1 lc rgb "#02D8DF" t '$\ev{K(\vb*{r}_i)}+\ev{V_{lj}(\vb*{r}_i)}$',\
    "../birthdeath/valencum.dat" u 0:($1*10.22/64) w l lw 1.5 lc rgb "#0000FF" t 'Cum. avg. of $\ev{K(\vb*{r}_i)}+\ev{V_{lj}(\vb*{r}_i)}$'