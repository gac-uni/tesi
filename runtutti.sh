#!/bin/bash
echo inizio esecuzione

echo inizio varmc
cd varmc
start=`date +%s`
./varmc_mt.out
end=`date +%s`
runtime_varmc=$((end-start))
echo tempo di esecuzione varmc $runtime_varmc
cd ../

echo inizio varmcnomrt
cd varmcnomrt
start=`date +%s`
./varmc_mt.out
end=`date +%s`
runtime_varmcnomrt=$((end-start))
echo tempo di esecuzione varmcnomrt $runtime_varmcnomrt
cd ../

echo inizio varmcnomrt_detbal
cd varmcnomrt_detbal
start=`date +%s`
./varmc_mt.out
end=`date +%s`
runtime_varmcnomrt_detbal=$((end-start))
echo tempo di esecuzione varmcnomrt_detbal $runtime_varmcnomrt_detbal
cd ../