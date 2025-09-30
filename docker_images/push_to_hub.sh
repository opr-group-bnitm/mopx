#!/usr/bin/env bash

images=(canu common medaka flye structural_variants)
docker_user=oprgroup

for img in "${images[@]}";
do
    version=$(cat $img/version.txt)
    echo "mopx_$img:$version"
    docker push "${docker_user}/mopx_${img}:$version"
    echo ""
done
