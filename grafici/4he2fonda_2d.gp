set terminal epslatex font "Times-New-Roman, 12"
#set terminal qt persist
set out "4he2fonda_2d.tex"

set isosamples 70
set surface

unset border
unset key 
unset grid

set tics nomirror

unset xtics
set xtics format ""
set xtics scale 0
unset ytics
set ytics format ""
set ytics scale 0
unset ztics
set ztics format ""
set ztics scale 0

set label "|" at -8,0,0 font "Times, 5"
set label "$-8$" at -8.5,0,-0.1 font "Times, 10"
# set label "|" at -6,0,0 font "Times, 5"
# set label "$-6$" at -6.5,0,-0.1 font "Times, 10"
set label "|" at -4,0,0 font "Times, 5"
set label "$-4$" at -4.5,0,-0.1 font "Times, 10"
# set label "|" at -2,0,0 font "Times, 5"
# set label "$-2$" at -2.5,0,-0.1 font "Times, 10"
set label "|" at 0,0,0 font "Times, 5"
set label "$0$" at -0.5,0,-0.1 font "Times, 10"
# set label "|" at 2,0,0 font "Times, 5"
# set label "$2$" at 1.5,0,-0.1 font "Times, 10"
set label "|" at 4,0,0 font "Times, 5"
set label "$4$" at 3.5,0,-0.1 font "Times, 10"
# set label "|" at 6,0,0 font "Times, 5"
# set label "$6$" at 5.5,0,-0.1 font "Times, 10"
set label "|" at 8,0,0 font "Times, 5"
set label "$8$" at 7.5,0,-0.1 font "Times, 10"

set label "|" at 0,-8,0 font "Times, 5"
set label "$-8$" at 0,-8.5,-0.1 font "Times, 10"
# set label "|" at 0,-6,0 font "Times, 5"
# set label "$-6$" at 0,-6.5,-0.1 font "Times, 10"
set label "|" at 0,-4,0 font "Times, 5"
set label "$-4$" at 0,-4.5,-0.1 font "Times, 10"
# set label "|" at 0,-2,0 font "Times, 5"
# set label "$-2$" at 0,-2.5,-0.1 font "Times, 10"
# set label "|" at 0,2,0 font "Times, 5"
# set label "$2$" at 0,1.5,-0.1 font "Times, 10"
set label "|" at 0,4,0 font "Times, 5"
set label "$4$" at 0,3.5,-0.1 font "Times, 10"
# set label "|" at 0,6,0 font "Times, 5"
# set label "$6$" at 0,5.5,-0.1 font "Times, 10"
set label "|" at 0,8,0 font "Times, 5"
set label "$8$" at 0,7.5,-0.1 font "Times, 10"

set label "-" at 0,0,0.5 font "Times, 10"
set label "$0.5$" at 0.5,0,0.5 font "Times, 10"
set label "-" at 0,0,1 font "Times, 10"
set label "$1$" at 0.5,0,1 font "Times, 10"
set label "-" at 0,0,1.5 font "Times, 10"
set label "$1.5$" at 0.5,0,1.5 font "Times, 10"

set label "$r_x [\\text{\\AA}]$" at 3, -10, 0
set label "$r_y [\\text{\\AA}]$" at -10, 2, 0
set label "$\\Psi(\\vb*{r})$" at 3, 3, 1.7 


set style arrow 1 head filled size char 1,20,50

set arrow 1 from -10,0,0 to 10,0,0 arrowstyle 1
set arrow 2 from 0,-10,0 to 0,10,0 arrowstyle 1
set arrow 3 from 0,0,-0.3 to 0,0,1.6 arrowstyle 1


f(x, y)=exp(-0.5*(2.5/(x**2+y**2)**(0.5))**5)

splot f(x, y) lc rgb "#0000FF"