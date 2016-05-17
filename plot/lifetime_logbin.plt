set terminal png
set output "lifetime_logbin.png"

unset key
set xlabel "Lifetime"
set ylabel "Frequency"

set format x "10^{%L}"
set logscale x
set format y "10^{%L}"
set logscale y

set style data lp

# stretched exponential fit
c = 2.5
t = 4.0
d = 0.5
h(x)=c+log10( exp(-(10.0**x/t)**d) )
i(x)=10.0**c*exp(-(x/t)**d)
fit [*:*] h(x) "lifetime_logbin.dat" u (log10($1)):(log10($2)) via c,t,d

p "lifetime_logbin.dat" u 1:2 title "Simulation", i(x) title "Stretched-exp. fit"

