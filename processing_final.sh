#!/bin/bash
cd processing_final
workdir="../"$1
outname="finalvalues_"$1".txt"
echo "analisi finale dei dati contenuti in "$1 > $outname
echo "" >> $outname 

# --- variabili che non do in input ---
pinitken=0
pinitpot=0
pinitstim1=0
pinitstim2=0
npassiwalker=500000
nwalkers=40
nnparticelle=64

# --- dimensionalizzo le energie ---
if [[ "$2" == "dim" ]];then
    echo dimensionalizzo
    ./dimensionalizza.sh $workdir $nparticelle
    echo fatto
fi

# --- calcolo le autocorrelazioni ---
echo autocorrelazioni
echo " --- output autocorr_mt.out --- " >> $outname
echo "" >> $outname
./autocorr_mt.out $workdir pinitken pinitpot pinitstim1 pinitstim2 | tee -a $outname 
echo fatto
sed -i '5d;6d;7d;8d' $outname
echo "" >> $outname
# --- prendo gli 1+t in input ---
echo 1+tau di ken
read v1ptken
echo 1+tau di pot
read v1ptpot
echo 1+tau di stim1
read v1ptstim1
echo 1+tau di stim2
read v1ptstim2
echo "valori di 1+tau dati in input--- ken: "$v1ptken", pot: "$v1ptpot", stim1: "$v1ptstim1", stim2: "$v1ptstim2 | tee -a $outname
echo "" >> $outname

# --- calcolo medie e varianze su tutti gli walkers --- 
echo medievarianze
echo " --- output medievarianze.out --- " >> $outname
echo "" >> $outname
./medievarianze.out $workdir $v1ptken $v1ptpot $v1ptstim1 $v1ptstim2 $nwalkers $npassiwalker | tee -a $outname
echo fatto




