
set xlabel 'k'
set ylabel 'Vijk'

set xzeroaxis

set title ' '

set key right center



plot 'k_vs_Vijk_tol1.dat' u 1:2 t 'tol=1E-4' with p  pointsize 1 pointtype 7 lc rgb '0x2c115f','k_vs_Vijk_tol2.dat' u 1:2 t 'tol=1E-5' with p  pointsize 1 pointtype 7 lc rgb '0x721f81', 'k_vs_Vijk_tol3.dat' u 1:2 t 'tol=1E-6' with p  pointsize 1 pointtype 7 lc rgb '0xb73779', 'k_vs_Vijk_tol4.dat' u 1:2 t 'tol=1E-7' with p  pointsize 1 pointtype 7 lc rgb '0xf1605d', 'k_vs_Vijk_tol5.dat' u 1:2 t 'tol=1E-8' with p  pointsize 1 pointtype 7 lc rgb '0xfeb078'

set terminal png
set output 'grafico_k_vs_Vijk.png'
replot
exit
