#!/usr/bin/env bash

images=(canu common medaka)

for img in "${images[@]}";
do
    echo $img
    docker build --no-cache -t mpox_analyse_${img}:0.0.1 $img/
    echo ""
done
