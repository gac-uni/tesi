#include <iostream>
#include <cassert>
#include "../readcsv.h"
#include "../gpd.h"

using namespace std;

void blocking(int piniziale, int npuntiplot, int npassi, double * vettoreiniziale, double * vectlarghezzablocchi, int * vectnblocchi, double * varblocchi, double * dFfinale){
    for(int i=0;i<npuntiplot;i++){
        int nblocchi=vectnblocchi[i];
        int larghezzablocco = vectlarghezzablocchi[i];
        double var=0;
        double df=0;
        double mediatot=0;
        double mediasb[nblocchi]={0};
        for(int j=0;j<nblocchi;j++){
            for(int k=0;k<larghezzablocco;k++){
                mediasb[j]+=vettoreiniziale[piniziale+larghezzablocco*j+k];
            }
            mediasb[j]=mediasb[j]/(double)larghezzablocco;
            mediatot+=mediasb[j];
        }
        mediatot=mediatot/(double)nblocchi;
        for(int j=0;j<nblocchi;j++){
            var+=(mediasb[j]-mediatot)*(mediasb[j]-mediatot);
        }
        var=var/(double)nblocchi;
        dFfinale[i]=sqrt(var/(double)nblocchi);
        varblocchi[i] = var;
    }
}

void askuserdf(double * vectdf, double * vectlblocchi, int npuntiplot, double * df, double * lblocco, string label){
    GnuplotDriver * gp;
    gp = new GnuplotDriver;
    gp->p();
    gp->t(label);
    gp->y("\\Delta F");
    gp->x("larghezza blocchi");
    gp->plot(vectlblocchi, vectdf, npuntiplot);
    delete gp; 
    cout << "df e larghezzablocco di"+label<<endl;
    cin >> *df >> *lblocco;
}

double * inputdata(string filename, int * npassi){
    FILE * fp = fopen(filename.c_str(), "r");
    readcsv pd(fp);
    double * v = new double[pd.nrighe];
    double * vtemp = new double[pd.nrighe];
    pd.dati(&vtemp);
    for(int i=0;i<pd.nrighe;i++){
        v[i]=vtemp[i];
    }
    if(*npassi == -1){
        *npassi = pd.nrighe;
    }else{
        assert(pd.nrighe == *npassi);
    }
    fclose(fp);
    return v;
}

void medvar(double * v, int npassi, double * media, double * var){
    *media=0;
    for(int i=0;i<npassi;i++){
        *media+=v[i];
    }
    *media/=(double)npassi;
    *var=0;
    for(int i=0;i<npassi;i++){
        *var+=(v[i]-*media)*(v[i]-*media);
    }
    *var/=(double)npassi;
}

int main(int argc, char ** argv){
    string workdir;
    int pinizialeken=0;
    int pinizialestim1=0;
    int pinizialestim2=0;
    int pinizialepot=0;
    switch(argc){
        case 2:
            workdir=argv[1];
            break;
        case 6:
            workdir=argv[1];
            sscanf(argv[2], "%d", &pinizialeken);
            sscanf(argv[3], "%d", &pinizialepot);
            sscanf(argv[4], "%d", &pinizialestim1);
            sscanf(argv[5], "%d", &pinizialestim2);

    }

    // --- input dei dati ---
    int npassiken=-1;
    int npassistim1=-1;
    int npassistim2=-1;
    int npassipot=-1;
    double * vKen=inputdata(workdir+"/mt_val_dim/ken_0.dat", &npassiken);
    double * vPot=inputdata(workdir+"/mt_val_dim/pot_0.dat", &npassistim1);
    double * vstim1=inputdata(workdir+"/mt_val_dim/stim1_0.dat", &npassistim2);
    double * vstim2=inputdata(workdir+"/mt_val_dim/stim2_0.dat", &npassipot); 

    npassiken -= pinizialeken; // npassi tenendo conto che scarto i primi piniziale per l'equilibrazione
    npassistim1 -= pinizialestim1;
    npassistim2 -= pinizialestim2;
    npassipot -= pinizialepot;

    int npuntiplot = 1000;
    int deltalarghezzaken=floor(npassiken/(5*npuntiplot+1))+0.5;
    int deltalarghezzapot=floor(npassipot/(5*npuntiplot+1))+0.5;
    int deltalarghezzastim1=floor(npassistim1/(5*npuntiplot+1))+0.5;
    int deltalarghezzastim2=floor(npassistim2/(5*npuntiplot+1))+0.5;

    double finken[npuntiplot];
    double varbken[npuntiplot];
    double finpot[npuntiplot];
    double varbpot[npuntiplot];
    double finstim1[npuntiplot];
    double varbstim1[npuntiplot];
    double finstim2[npuntiplot];
    double varbstim2[npuntiplot];
    int nblocchiken[npuntiplot];
    double larghezzablocchiken[npuntiplot];
    int nblocchistim1[npuntiplot];
    double larghezzablocchistim1[npuntiplot];
    int nblocchistim2[npuntiplot];
    double larghezzablocchistim2[npuntiplot];
    int nblocchipot[npuntiplot];
    double larghezzablocchipot[npuntiplot];
    for(int i=0;i<npuntiplot;i++){
        larghezzablocchiken[i] = deltalarghezzaken*(1+i);
        nblocchiken[i]=floor(((double)npassiken)/(double)larghezzablocchiken[i])+0.5;
        larghezzablocchipot[i] = deltalarghezzapot*(1+i);
        nblocchipot[i]=floor(((double)npassipot)/(double)larghezzablocchipot[i])+0.5;
        larghezzablocchistim1[i] = deltalarghezzastim1*(1+i);
        nblocchistim1[i]=floor(((double)npassistim1)/(double)larghezzablocchistim1[i])+0.5;
        larghezzablocchistim2[i] = deltalarghezzastim2*(1+i);
        nblocchistim2[i]=floor(((double)npassistim2)/(double)larghezzablocchistim2[i])+0.5;
    }

    blocking(pinizialeken, npuntiplot, npassiken, vKen, larghezzablocchiken, nblocchiken, varbken, finken);
    blocking(pinizialepot, npuntiplot, npassipot, vPot, larghezzablocchipot, nblocchipot, varbpot, finpot);
    blocking(pinizialestim1, npuntiplot, npassistim1, vstim1, larghezzablocchistim1, nblocchistim1, varbstim1, finstim1);
    blocking(pinizialestim2, npuntiplot, npassistim2, vstim2, larghezzablocchistim2, nblocchistim2, varbstim2, finstim2);

    double lbloccoken;
    double dfken;
    double lbloccopot;
    double dfpot;
    double lbloccostim1;
    double dfstim1;
    double lbloccostim2;
    double dfstim2;

    askuserdf(finken, larghezzablocchiken, npuntiplot, &dfken, &lbloccoken, "ken");
    askuserdf(finpot, larghezzablocchipot, npuntiplot, &dfpot, &lbloccopot, "pot");
    askuserdf(finstim1, larghezzablocchistim1, npuntiplot, &dfstim1, &lbloccostim1, "stim1");
    askuserdf(finstim2, larghezzablocchistim2, npuntiplot, &dfstim2, &lbloccostim2, "stim2");

    // calcolo delle varianze vere e delle medie
    double mediaken;
    double varken;
    double mediapot;
    double varpot;
    double mediastim1;
    double varstim1;
    double mediastim2;
    double varstim2;
    medvar(vKen, npassiken, &mediaken, &varken);
    medvar(vPot, npassipot, &mediapot, &varpot);
    medvar(vstim1, npassistim1, &mediastim1, &varstim1);
    medvar(vstim2, npassistim2, &mediastim2, &varstim2);

    // FILE * out; 
    // string nome =workdir+"/mt_val_dim/varianze_autocorr.dat"; 
    // out = fopen(nome.c_str(), "w");
    printf("ken --- media: %lf, errore: %lf, var_corretta: %lf\n", mediaken, sqrt(varken/(double)npassiken), dfken);
    printf("pot --- media: %lf, errore: %lf, var_corretta: %lf\n", mediapot, sqrt(varpot/(double)npassipot), dfpot);
    printf("stim1 --- media: %lf, errore: %lf, var_corretta: %lf\n", mediastim1, sqrt(varstim1/(double)npassistim1), dfstim1);
    printf("stim2 --- media: %lf, errore: %lf, var_corretta: %lf\n", mediastim2, sqrt(varstim2/(double)npassistim2), dfstim2);
    // fprintf(out, "ken --- media: %lf, errore: %lf, var_corretta: %lf\n", mediaken, sqrt(varken/(double)npassiken), dfken);
    // fprintf(out, "pot --- media: %lf, errore: %lf, var_corretta: %lf\n", mediapot, sqrt(varpot/(double)npassipot), dfpot);
    // fprintf(out, "stim1 --- media: %lf, errore: %lf, var_corretta: %lf\n", mediastim1, sqrt(varstim1/(double)npassistim1), dfstim1);
    // fprintf(out, "stim2 --- media: %lf, errore: %lf, var_corretta: %lf\n", mediastim2, sqrt(varstim2/(double)npassistim2), dfstim2);
    // fprintf(out, "\n");
    // calcolo stime del numero di blocchi usato per trovare 1+tau
    double stimanblocchiken = floor(npassiken/lbloccoken);
    double stimanblocchipot = floor(npassipot/lbloccopot);
    double stimanblocchistim1 = floor(npassistim1/lbloccostim1);
    double stimanblocchistim2 = floor(npassistim2/lbloccostim2);
    double sigquadfbstimaken = dfken*dfken*stimanblocchiken;
    double sigquadfbstimapot = dfpot*dfpot*stimanblocchipot;
    double sigquadfbstimastim1 = dfstim1*dfstim1*stimanblocchistim1;
    double sigquadfbstimastim2 = dfstim2*dfstim2*stimanblocchistim2;
    printf("ken --- (1+tau) = %lf, stima di df: %lf\n", (lbloccoken*sigquadfbstimaken/varken), sqrt(varken/npassiken*(lbloccoken*sigquadfbstimaken/varken)));
    printf("pot --- (1+tau) = %lf, stima di df: %lf\n", (lbloccopot*sigquadfbstimapot/varpot), sqrt(varpot/npassipot*(lbloccopot*sigquadfbstimapot/varpot)));
    printf("stim1 --- (1+tau) = %lf, stima di df: %lf\n", (lbloccostim1*sigquadfbstimastim1/varstim1), sqrt(varstim1/npassistim1*(lbloccostim1*sigquadfbstimastim1/varstim1)));
    printf("stim2 --- (1+tau) = %lf, stima di df: %lf\n", (lbloccostim2*sigquadfbstimastim2/varstim2), sqrt(varstim2/npassistim2*(lbloccostim2*sigquadfbstimastim2/varstim2)));
    // fprintf(out, "ken --- (1+tau) = %lf, stima di df: %lf\n", (lbloccoken*sigquadfbstimaken/varken), sqrt(varken/npassiken*(lbloccoken*sigquadfbstimaken/varken)));
    // fprintf(out, "pot --- (1+tau) = %lf, stima di df: %lf\n", (lbloccopot*sigquadfbstimapot/varpot), sqrt(varpot/npassipot*(lbloccopot*sigquadfbstimapot/varpot)));
    // fprintf(out, "stim1 --- (1+tau) = %lf, stima di df: %lf\n", (lbloccostim1*sigquadfbstimastim1/varstim1), sqrt(varstim1/npassistim1*(lbloccostim1*sigquadfbstimastim1/varstim1)));
    // fprintf(out, "stim2 --- (1+tau) = %lf, stima di df: %lf\n", (lbloccostim2*sigquadfbstimastim2/varstim2), sqrt(varstim2/npassistim2*(lbloccostim2*sigquadfbstimastim2/varstim2)));

    // fclose(out);


    return 0;

}