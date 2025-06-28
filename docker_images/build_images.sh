#!/usr/bin/env bash

images=(canu common medaka flye structural_variants)

for img in "${images[@]}";
do
    echo $img
    docker build --no-cache -t mopx_${img}:0.0.1 $img/
    echo ""
done
