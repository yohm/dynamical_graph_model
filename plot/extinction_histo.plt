set terminal png
set output "extinction_histo.png"

unset key
set xlabel "Extinction Size"
set ylabel "Frequency"

set format y "10^{%L}"
set logscale y

set style data lp
p "extinction_histo.dat" u 1:2

