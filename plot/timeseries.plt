set terminal png
set output "timeseries.png"

unset key
set xlabel "Time (k steps)"
set ylabel "Number of Species"
set style data lp
p "timeseries.dat" u ($1/1000):2

