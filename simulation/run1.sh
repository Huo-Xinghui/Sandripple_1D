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
#run_simulation "uStar030_240log50_0_2650_300" "1.out"
#run_simulation "uStar040_240log50_0_2650_300" "2.out"
#run_simulation "uStar050_240log50_0_2650_300" "3.out"
run_simulation "uStar035_269log100_0_2650_300" "4.out"
run_simulation "uStar045_269log100_0_2650_300" "5.out"
run_simulation "uStar055_269log100_0_2650_300" "6.out"
cd "/home/ekalhxh/2025/Q_on_flat_bed"
