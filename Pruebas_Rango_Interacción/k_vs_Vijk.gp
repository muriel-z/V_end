
set xlabel 'k'
set ylabel 'Vijk'

set xzeroaxis

set title ' '

set key right center



plot 'k_vs_Vijk.dat' u 1:2 t 'V' with p  pointsize 1 pointtype 7 lc rgb '0x2c115f'

set terminal png
set output 'grafico_k_vs_Vijk.png'
replot
exit
