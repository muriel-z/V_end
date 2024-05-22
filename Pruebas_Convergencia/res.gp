
set xlabel '# de iteracion'
set ylabel 'Res (V0(i,j,k)-V(i,j,k))**2/(nvx*nvy*nvz) '

set xzeroaxis
set logscale y

set title ' '

set key off


plot 'res_vs_it.dat' u 1:2 with p  pointsize 1 pointtype 7

set terminal png
set output 'grafico_res_vs_it.png'
replot
exit

