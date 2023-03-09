#include "variables.h"

// --- file ---
std::string valdir="varmc_val/";
std::string plotdir = "varmc_plots/";

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
int npassi=10000;
int freqcampionamento=1;
double tau=0.0001;
double a1=2.5/sigma;
double a2=5;

// --- vettori costruiti con i parametri del sistema ---
double parametripot[2]={1, 1};
double parametrilogprob[2]={a1, a2};

// --- parametri calcolati dai parametri del sistema ---
int nce=(nparticelle/particelleperce);
double volce=(double)nparticelle/(densita*(double)nce);
double latoce=pow(volce, 1./3.);
double latoscatola=latoce*pow(nce, 1./3.);
int npassiplot = npassi/freqcampionamento+0.5; 
    