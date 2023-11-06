
set xlabel 'malla en z'
set ylabel 'Res sum '

set xzeroaxis
set logscale y

set title ' '

set key off


plot 'res_sum_100.dat' u 3:4:0 with p  pointsize 1.5 pointtype 7

set terminal png
set output 'grafico_res_sum_100_z.png'
replot
exit

