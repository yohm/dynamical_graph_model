#!/bin/bash -ex

script_dir=$(cd $(dirname $BASH_SOURCE); pwd)
$script_dir/main.out $@

for pyfile in $( ls $script_dir/plot/*.py ); do
  python "$pyfile"
done

