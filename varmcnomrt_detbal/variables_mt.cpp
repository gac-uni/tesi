#include "variables_mt.h"

// --- file ---
FILE * out = fopen("dmc_out/out.dat", "w");
FILE * valori=fopen("dmc_out/valori.dat", "w");

// --- parametri per i reticoli della posizione iniziale ---
// int particelleperce=4; // particelle per cella elementare
// double posizioniince[4][3]={{0, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}}; // posizione della particella nella cella elementare
int particelleperce=1; // particelle per cella elementare
double posizioniince[1][3]={{0,0,0}}; // posizione della particella nella cella elementare

// --- parametri del sistema ---
int nparticelle=64;
double eps=10.22;
double sigma=2.556;
double densita=0.0218*sigma*sigma*sigma; // in unità adimensionali \tilda(r)=r/sigma
double a1=2.5/sigma;
double a2=5;

// --- parametri variazionali ---
double tau=0.00005;

// --- parametri defusion ---
int npassiblock = 500000;
int freqcampionamento=1;
int nwalkers = 40;
// int nblocks = 1;

// --- vettori costruiti con i parametri del sistema ---
double parametripot[2]={1, 1};
double parametrilogprob[2]={a1, a2};

// --- parametri calcolati dai parametri del sistema ---
int nce=(nparticelle/particelleperce);
double volce=(double)nparticelle/(densita*(double)nce);
double latoce=pow(volce, 1./3.);
double latoscatola=latoce*pow(nce, 1./3.);
int npassiplot = npassiblock/freqcampionamento+0.5; 

// --- variabili multithreading ---
int n_threads=std::thread::hardware_concurrency()-1;
    