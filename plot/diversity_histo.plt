set terminal png
set output "diversity_histo.png"

unset key
set xlabel "Number of Species"
set ylabel "Frequency"

set format y "10^{%L}"
set logscale y

set style data lp
p "diversity_histo.dat" u 1:2

