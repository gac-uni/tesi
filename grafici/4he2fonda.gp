set terminal epslatex font "Times-New-Roman, 12"
# set terminal qt persist
set out "4he2fondapresentazione.tex"
set samples 1000

unset border
set key at 17, 0.4
unset grid

set xtics axis
set ytics axis

set style arrow 1 head filled size char 1,20,50

set arrow 1 from -0.5,0 to 13,0 arrowstyle 1
set arrow 2 from 0,-0.08 to 0, 1.2 arrowstyle 1

set xrange[0:13]
set yrange[0:1]

set tmargin 6
set rmargin 10

set tics scale 0.5

set xlabel '$\qty|\vb*{r}| \quad [\unit{\angstrom}]$'
set ylabel '$\Psi(\vb*{r})$' 

# set label '$\qty|\vb*{r}|$' at 12, -0.2
# set label '$\Psi(\vb*{r})$' at -2.5, 1.1 



f(x)=exp(-0.5*(2.5/(abs(x))**(0.5))**5)

plot f(x) lc rgb "#0000FF" t '$ \Psi_{\alpha}(\vb*{r}) = \exp[-\frac{1}{2}\qty(\frac{\alpha_1}{r})^{\alpha_2}]$'