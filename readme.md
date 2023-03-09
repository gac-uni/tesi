## Metodi VMC
nelle cartelle varmc, varmcnomrt, varmcnomrt_detbal si trovano rispettivamente i codici per il metodo variazionale con metropolis, quello con la drift velocity e quello con la drift velocity e la correzione con accept/reject

in ciascuna di queste cartella varmc.cpp è il sorgente per l'evoluzione di un walker singolo mentre varmc_multithread.cpp è il sorgente per l'evoluzione di nwalkers walkers, parallelizzato.
varmc_val contiene i risultati di varmc.cpp e mt_val contiene i risultati di varmc_multithread.cpp.
processing_final.sh calcola i valori finali sui risultati di varmc_multithread.cpp. La cartella a cui si riferisce (varmc, varmcnomrt, varmcnomrt_detbal) va esplicitata come argomento di processing_final.sh. Per far funzionare processing_final.sh vanno compilati i file della cartella processing_final: autocorr_mt.cpp -> autocorr_mt.out, medievarianze.cpp -> medievarianze.out.

## Metodi DMC

### dmc_postp
getvarianze.cpp calcola le varianze corrette con le autocorrelazioni (trovate da processing_final.sh) e le mette in un file. dmc_postp.cpp richiede queste varianze come array (gia ci sono quelli usati per la tesi) e produce i risultati per dmc_postp.

### birthdeath

birthdeath.cpp esegue l'algoritmo di birth-death.