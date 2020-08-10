#!/bin/bash
set -ex 

matlab -nodisplay -r \
    "addpath(genpath('./data')); Arch_Efficiency_Simulations" &

for dataset in simulation_workspace_M{27,32,38}; do

matlab -nodisplay -r \
    "addpath('./data'); FalseNegPosRate_Simulations('$dataset')" &

done
wait