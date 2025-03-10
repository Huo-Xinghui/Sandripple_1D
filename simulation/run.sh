#!/bin/bash

file="250225.f90"
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

    cd "/home/ekalhxh/ripple/coll11/$folder"
    rm -f nohup.out
    # uStar
    line_number=59
    sed -i "${line_number}s/[0-9]\+\.[0-9]\{2\}/$uStar/" "$file"
    # dpa
    line_number=33
    sed -i "${line_number}s/[0-9]\+\.[0-9]\+e-[0-9]\+/$dpa/" "$file"
    # rho
    line_number=60
    sed -i "${line_number}s/[0-9]\+\.[0-9]\{3\}/$rho/" "$file"
    # rhoP
    line_number=43
    sed -i "${line_number}s/[0-9]\{4\}\.[0-9]/$rhoP/" "$file"
    # endTime
    line_number=77
    sed -i "${line_number}s/[0-9]\+.[0-9]/$endTime/" "$file"
    # nNodes
    line_number=17
    sed -i "${line_number}s/[0-9]\+/$nNodes/" "$file"
    # stdd
    line_number=34
    sed -i "${line_number}s/[0-9]\+\.[0-9]\+e-[0-9]\+/$stdd/" "$file"

    mpif90 -o "$output" "$file"
    nohup mpirun --bind-to none -np $nNodes "$output" &
}
run_simulation "uStar050_200_0_2650_600" 0.50 2.0e-4 1.263 2650.0 600.0 1.0e-4 "200.out"
run_simulation "uStar050_400_0_2650_600" 0.50 4.0e-4 1.263 2650.0 600.0 1.0e-4 "400.out"
run_simulation "uStar050_500_0_2650_600" 0.50 5.0e-4 1.263 2650.0 600.0 1.0e-4 "500.out"
run_simulation "uStar050_300_2_2650_600" 0.50 3.0e-4 2.0 2650.0 600.0 1.0e-4 "2.out"
run_simulation "uStar050_300_0_1000_600" 0.50 3.0e-4 1.263 1000.0 600.0 1.0e-4 "1000.out"
run_simulation "uStar050_300_0_4000_600" 0.50 3.0e-4 1.263 4000.0 600.0 1.0e-4 "4000.out"
cd "/home/ekalhxh/ripple/coll11"
