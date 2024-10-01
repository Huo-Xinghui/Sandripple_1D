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
    local output=$7

    cd "/home/ekalhxh/sandripple/ripple/coll/$folder"
    rm -f nohup.out
    # uStar
    line_number=56
    sed -i "${line_number}s/[0-9]\+\.[0-9]\{2\}/$uStar/" "$file"
    # dpa
    line_number=32
    sed -i "${line_number}s/[0-9]\+\.[0-9]\+e-[0-9]\+/$dpa/" "$file"
    # rho
    line_number=57
    sed -i "${line_number}s/[0-9]\+\.[0-9]\{3\}/$rho/" "$file"
    # rhoP
    line_number=40
    sed -i "${line_number}s/[0-9]\{4\}\.[0-9]/$rhoP/" "$file"
    # endTime
    line_number=73
    sed -i "${line_number}s/[0-9]\+.[0-9]/$endTime/" "$file"
    # nNodes
    line_number=17
    sed -i "${line_number}s/[0-9]\+/$nNodes/" "$file"

    mpif90 -o "$output" "$file"
    nohup mpirun --bind-to none -np $nNodes "$output" &
}
run_simulation "uStar040_150and350_0_2650_3600" 0.40 2.5e-4 1.263 2650.0 3600.0 "040.out"
run_simulation "uStar045_150and350_0_2650_3600" 0.45 2.5e-4 1.263 2650.0 3600.0 "045.out"
run_simulation "uStar050_150and350_0_2650_3600" 0.50 2.5e-4 1.263 2650.0 3600.0 "050.out"
run_simulation "uStar055_150and350_0_2650_3600" 0.55 2.5e-4 1.263 2650.0 3600.0 "055.out"
run_simulation "uStar060_150and350_0_2650_3600" 0.60 2.5e-4 1.263 2650.0 3600.0 "060.out"
cd "/home/ekalhxh/sandripple/ripple/coll"