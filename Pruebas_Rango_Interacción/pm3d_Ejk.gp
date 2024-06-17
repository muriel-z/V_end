set xlabel ''
set ylabel ''

set xzeroaxis


set title ' '

set key bottom left

set pm3d map


splot 'Ejk.dat' u 1:2:3

set terminal png 
set output 'grafico_Ejk.png'
replot 
exit
