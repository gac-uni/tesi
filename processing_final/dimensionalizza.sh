#!/bin/bash

workdir=$1
nparticelle=$2
param_energia=10.22
cost_dim=$(bc <<< "scale=3; $param_energia/64")
mkdir -p $workdir"/mt_val_dim"
for i in {0..39}
do
    fken="ken_"$i".dat"
    fpot="pot_"$i".dat"
    fstim1="stim1_"$i".dat"
    fstim2="stim2_"$i".dat"
    awk -vcost_dim=$cost_dim '{temp = $0*cost_dim; print temp}' $workdir"/mt_val/"$fken > $workdir"/mt_val_dim/"$fken
    awk -vcost_dim=$cost_dim '{temp = $0*cost_dim; print temp}' $workdir"/mt_val/"$fpot > $workdir"/mt_val_dim/"$fpot
    awk -vcost_dim=$cost_dim '{temp = $0*cost_dim; print temp}' $workdir"/mt_val/"$fstim1 > $workdir"/mt_val_dim/"$fstim1
    awk -vcost_dim=$cost_dim '{temp = $0*cost_dim; print temp}' $workdir"/mt_val/"$fstim2 > $workdir"/mt_val_dim/"$fstim2 
done