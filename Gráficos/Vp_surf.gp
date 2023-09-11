
set xlabel 'malla en z'
set ylabel 'malla en y '

set xzeroaxis


set title 'Valores del voltaje'

set key off


splot 'Vp.dat' u 4:3:5 with points palette pointsize 1.5 pointtype 7

set terminal png
set output 'grafico_Vp_surf.png'
replot 
exit

