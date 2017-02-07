#!/bin/bash

set -ex

script_dir=$(cd $(dirname $0); pwd)
cd "$script_dir"
make

if [ -n "$OACIS_ROOT" ];
then
  "$OACIS_ROOT/bin/oacis_ruby" "$script_dir/register_oacis.rb"
fi
