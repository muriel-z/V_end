
set xlabel 'k'
set ylabel 'Res (V0(i,j,k)-V(i,j,k))**2 '

set xzeroaxis
set logscale y

set title ' '

set key off


plot 'res_100.dat' u 3:4:0 with p  pointsize 1.5 pointtype 7

set terminal png
set output 'grafico_res_z_100.png'
replot
exit

