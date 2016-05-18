set terminal png
set output "timeseries.png"

set key
set xlabel "Time (k steps)"
set ylabel "Number of Species"
set y2tics
set y2label "Clustering Coefficient"
set style data lp
p "timeseries.dat" u ($1/1000):2 title "# of Species", "" u ($1/1000):3 axis x1y2 title "CC"

