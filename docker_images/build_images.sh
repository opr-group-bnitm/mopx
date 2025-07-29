#!/usr/bin/env bash

images=(canu common medaka flye structural_variants)

for img in "${images[@]}";
do
    version=$(cat $img/version.txt)
    echo "mopx_$img:$version"
    docker build --no-cache -t "oprgroup/mopx_${img}:$version" $img/
    echo ""
done
