#!/usr/bin/gnuplot

set xlabel 'malla en z'
set ylabel 'malla en y '

set xzeroaxis


set title 'Valores del voltaje'

set key off

set hidden3d 

splot [][][290:300] 'Vp.dat' u 2:3:4:5 with points linecolor palette pointsize 1.5 pointtype 7

pause -1

      /* "<sed '1,2d' pos_final.xyz"  with points lc rgb "red" */

/* set terminal png */
/* set output 'grafico_Vp_surf.png' */
/* replot  */
/* exit */

