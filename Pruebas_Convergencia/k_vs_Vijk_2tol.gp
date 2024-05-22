
set xlabel 'k'
set ylabel 'Vijk'

set xzeroaxis

set title ' '

set key right center


plot 'k_vs_Vijk_tol5.dat' u 1:2 t 'tol=1E-8' with p  pointsize 1 pointtype 7 lc rgb '0x2c115f','k_vs_Vijk_tol6.dat' u 1:2 t 'tol=1E-9' with p  pointsize 1 pointtype 7 lc rgb '0x721f81'

set terminal png
set output 'grafico_k_vs_Vijk_tol5y6.png'
replot
exit