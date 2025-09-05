#!/bin/bash

file="250814.f90"
run_simulation() {
    local folder=$1
    local output=$2

    cd "/home/ekalhxh/2025/Q_on_flat_bed/$folder"
    rm -f nohup.out
    mpif90 -o "$output" "$file"
    nohup mpirun --bind-to none -np 5 "$output" &
}
run_simulation "uStar050_290log97_0_2650_300" "3.out"
run_simulation "uStar060_290log97_0_2650_300" "4.out"
run_simulation "uStar030_197log65_0_2650_300" "5.out"
run_simulation "uStar040_197log65_0_2650_300" "1.out"
run_simulation "uStar050_197log65_0_2650_300" "2.out"
run_simulation "uStar060_197log65_0_2650_300" "6.out"
cd "/home/ekalhxh/2025/Q_on_flat_bed"
