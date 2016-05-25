
set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0

set multiplot
set size 0.5,0.5

set origin 0.0,0.5
set grid
set key
set title "Temperature along the axis"
set xlabel "axial coordinate [mm]"
set ylabel "temperature [K]"
set xrange [0:140]
set yrange [200:1700]

plot \
  '02-steady-state/postProcessing/sets/4000/axis_p_T_H2_O2_H2O_N2.xy' u ($1*1000):3 w l t 'coarse', \
  '03-fine-grid/postProcessing/sets/10000/axis_p_T_H2_O2_H2O_N2.xy' u ($1*1000):3 w l t 'fine', \
  'exp/axis-T.exp' u 1:2 w p t 'Exp. (CARS)', \
  'exp/axis-T.exp' u 3:4 w p t 'Exp. (Raman)'

set origin 0.5,0.0
set grid
set key
set title "Temperature along the radial coordinate @ 30 mm"
set xlabel "radial coordinate [mm]"
set ylabel "temperature [K]"
set xrange [0:20]
set yrange [200:2000]

plot \
  '02-steady-state/postProcessing/sets/4000/radial_30mm_p_T_H2_O2_H2O_N2.xy' u ($1*1000):3 w l t 'coarse', \
  '03-fine-grid/postProcessing/sets/10000/radial_30mm_p_T_H2_O2_H2O_N2.xy' u ($1*1000):3 w l t 'fine', \
  'exp/radial-30mm-T.exp' u 1:2 w p t 'Exp. (CARS)', \
  'exp/radial-30mm-T.exp' u 3:4 w p t 'Exp. (Raman)'

set origin 0.,0.
set grid
set key
set title "Temperature along the radial coordinate @ 20 mm"
set xlabel "radial coordinate [mm]"
set ylabel "temperature [K]"
set xrange [0:20]
set yrange [200:2000]

plot \
  '02-steady-state/postProcessing/sets/4000/radial_20mm_p_T_H2_O2_H2O_N2.xy' u ($1*1000):3 w l t 'coarse', \
  '03-fine-grid/postProcessing/sets/10000/radial_20mm_p_T_H2_O2_H2O_N2.xy' u ($1*1000):3 w l t 'fine', \
  'exp/radial-20mm-T.exp' u 1:2 w p t 'Exp. (CARS)', \
  'exp/radial-20mm-T.exp' u 3:4 w p t 'Exp. (Raman)'

set origin 0.5,0.5
set grid
set key
set title "Temperature along the radial coordinate @ 10 mm"
set xlabel "radial coordinate [mm]"
set ylabel "temperature [K]"
set xrange [0:20]
set yrange [200:2000]

plot \
  '02-steady-state/postProcessing/sets/4000/radial_10mm_p_T_H2_O2_H2O_N2.xy' u ($1*1000):3 w l t 'coarse', \
  '03-fine-grid/postProcessing/sets/10000/radial_10mm_p_T_H2_O2_H2O_N2.xy' u ($1*1000):3 w l t 'fine', \
  'exp/radial-10mm-T.exp' u 1:2 w p t 'Exp. (CARS)', \
  'exp/radial-10mm-T.exp' u 3:4 w p t 'Exp. (Raman)'

unset multiplot
