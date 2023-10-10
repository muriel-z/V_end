#!/usr/bin/gnuplot

set xlabel 'malla en z'
set ylabel 'malla en y '

set xzeroaxis


set title 'Valores del voltaje'

set key off

# set hidden3d 

# splot [][][] 'V.dat' u 1:2:3:4 with points linecolor palette pointsize 1.5 pointtype 7

splot [][][] "< awk  '($2==25)' V.dat" u 1:3:4 

pause -1

      /* "<sed '1,2d' pos_final.xyz"  with points lc rgb "red" */

/* set terminal png */
/* set output 'grafico_Vp_surf.png' */
/* replot  */
/* exit */

