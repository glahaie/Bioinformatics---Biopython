set xrange [-25:375]
unset key

plot 'contigs_taille.dat' u 1:2
min_y = GPVAL_DATA_Y_MIN
max_y = GPVAL_DATA_Y_MAX

f(x) = mean_y
fit f(x) 'contigs_taille.dat' u 1:2 via mean_y

# Plotting the minimum and maximum ranges with a shaded background
set label 1 gprintf("Minimum = %g", min_y) at 2, min_y-10
set label 2 gprintf("Maximum = %g", max_y) at 2, max_y+10
set label 3 gprintf("Moyenne = %g", mean_y) at 2, max_y+20
set label 4 "Taille des contigs" at 275, max_y+20
set term postscript
set output "contigs_taille.eps"
plot min_y with filledcurves y1=mean_y lt 1 lc rgb "#bbbbdd", \
max_y with filledcurves y1=mean_y lt 1 lc rgb "#bbddbb", \
'contigs_taille.dat' u 1:2 w p pt 2 lt 1 ps 1

