#!/usr/bin/env -S gnuplot -persist
set style data linespoints
# set style data lines
set key outside top center horizontal
set offsets graph 0.1, graph 0.1, graph 0.1, graph 0.1
set title "Instantaneous Profile"
set grid xtics
set grid ytics
set grid

process = "< awk '(NR>2){print;}' "
file = filename
plt_data = process.file
plot plt_data using 1:2 title "u_h", \
    # sin(pi * (x - 2.) - sin(pi * (x - 2.)) / pi) title "u_a"
	# exp(-8. * (x-1.) * (x-1.)) title "u_a"
	# "< awk '(NR>2){print;}' res.dat" using 1:4 title "p", \
	# "< awk '(NR>2){print;}' res.dat" using 1:2 title "ρ", \
	# "< awk '(NR>2){print;}' res.dat" using 1:3 title "u", \
	# "< awk '(NR>2){print;}' res.dat" using 1:5 title "e", \
	# "< awk '(NR>2){print;}' amp_0.01wn_0.1p.continue3000.txt" using 1:2 title "u"
pause mouse close
