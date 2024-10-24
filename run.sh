#!/bin/bash

file="240929.f90"
nNodes=5
run_simulation() {
    local folder=$1
    local uStar=$2
    local dpa=$3
    local rho=$4
    local rhoP=$5
    local endTime=$6
    local stdd=$7
    local output=$8

    cd "/home/ekalhxh/ripple/coll/$folder"
    rm -f nohup.out
    # uStar
    line_number=57
    sed -i "${line_number}s/[0-9]\+\.[0-9]\{2\}/$uStar/" "$file"
    # dpa
    line_number=32
    sed -i "${line_number}s/[0-9]\+\.[0-9]\+e-[0-9]\+/$dpa/" "$file"
    # rho
    line_number=58
    sed -i "${line_number}s/[0-9]\+\.[0-9]\{3\}/$rho/" "$file"
    # rhoP
    line_number=41
    sed -i "${line_number}s/[0-9]\{4\}\.[0-9]/$rhoP/" "$file"
    # endTime
    line_number=74
    sed -i "${line_number}s/[0-9]\+.[0-9]/$endTime/" "$file"
    # nNodes
    line_number=17
    sed -i "${line_number}s/[0-9]\+/$nNodes/" "$file"
    # stdd
    line_number=33
    sed -i "${line_number}s/[0-9]\+\.[0-9]\+e-[0-9]\+/$stdd/" "$file"


    mpif90 -o "$output" "$file"
    nohup mpirun --bind-to none -np $nNodes "$output" &
}
run_simulation "uStar040_300stdd100_0_2650_3600" 0.40 3.0e-4 1.263 2650.0 3600.0 1.0e-4 "40.out"
run_simulation "uStar045_300stdd100_0_2650_3600" 0.45 3.0e-4 1.263 2650.0 3600.0 1.0e-4 "45.out"
run_simulation "uStar055_300stdd100_0_2650_3600" 0.55 3.0e-4 1.263 2650.0 3600.0 1.0e-4 "55.out"
run_simulation "uStar060_300stdd100_0_2650_3600" 0.60 3.0e-4 1.263 2650.0 3600.0 1.0e-4 "60.out"
#run_simulation "uStar050_300stdd100_0_2650_3600" 0.50 3.0e-4 1.263 2650.0 3600.0 5.0e-6 "5.out"
cd "/home/ekalhxh/ripple/coll"