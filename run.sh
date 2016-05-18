#!/bin/bash -ex

script_dir=$(cd $(dirname $BASH_SOURCE); pwd)
$script_dir/main.out $1 $2 $3 $4

ruby $script_dir/plot/log_binning.rb lifetime.dat > lifetime_logbin.dat
for pltfile in $( ls $script_dir/plot/*.plt ); do
  gnuplot "$pltfile"
done

