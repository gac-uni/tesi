#include <iostream>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <thread>
#include "gpd.h"

#define NPASSIVAR 1000
#define NPUNTIPLOT rint()

class Walker{
    public:
        double peso;
        double * ken(){
            /*
            parametri: 
            0) a1
            1) a2
            */
            double * out;
            out = new double[3];
            /*
            out:
            0) ken
            1) stimatore1
            2) stimatore2
            */
            out[0]=0;
            out[1]=0;
            out[2]=0;
            for(int i=0;i<nparticelle; i++){
                double lap=0;
                double gradx=0;
                double grady=0;
                double gradz=0;
                double gradquad=0;
                for(int j=0;j<nparticelle;j++){
                    double dist=sqrt((ppvx[nparticelle*i+j]*ppvx[nparticelle*i+j])+(ppvy[nparticelle*i+j]*ppvy[nparticelle*i+j])+(ppvz[nparticelle*i+j]*ppvz[nparticelle*i+j]));
                    if(i!=j){
                        double laptemp=0;
                        laptemp+= -0.5*((ppvx[nparticelle*i+j]*ppvx[nparticelle*i+j])/(dist*dist))*(parametriken[0], parametriken[1])*parametriken[1]*(parametriken[1]+1)*(pow(dist, -parametriken[1]-2)+pow(lscatola-dist, -parametriken[1]-2));
                        laptemp+= -0.5*((ppvy[nparticelle*i+j]*ppvy[nparticelle*i+j])/(dist*dist))*(parametriken[0], parametriken[1])*parametriken[1]*(parametriken[1]+1)*(pow(dist, -parametriken[1]-2)+pow(lscatola-dist, -parametriken[1]-2));
                        laptemp+= -0.5*((ppvz[nparticelle*i+j]*ppvz[nparticelle*i+j])/(dist*dist))*(parametriken[0], parametriken[1])*parametriken[1]*(parametriken[1]+1)*(pow(dist, -parametriken[1]-2)+pow(lscatola-dist, -parametriken[1]-2));
                        if(j<i){
                            laptemp*=2;
                            gradx+= -0.5*(ppvx[nparticelle*i+j]/dist)*pow(parametriken[0], parametriken[1])*parametriken[1]*(-pow(dist, -parametriken[1]-1)+pow(lscatola-dist, -parametriken[1]-1));
                            grady+= -0.5*(ppvy[nparticelle*i+j]/dist)*pow(parametriken[0], parametriken[1])*parametriken[1]*(-pow(dist, -parametriken[1]-1)+pow(lscatola-dist, -parametriken[1]-1));
                            gradz+= -0.5*(ppvz[nparticelle*i+j]/dist)*pow(parametriken[0], parametriken[1])*parametriken[1]*(-pow(dist, -parametriken[1]-1)+pow(lscatola-dist, -parametriken[1]-1));
                        }                   
                        lap+=laptemp;        
                    }
                }
                gradquad=gradx*gradx+grady*grady+gradz*gradz;
                out[0]+= -6.0596*(lap+gradquad);
                out[1]+= -3.0298*lap;
                out[2]+= 6.0596*gradquad;
            }
            return out;
        }
        double pot(){
            /*
            parametri:
            0) eps
            1) sigma
            */
            double pot1=0;
            for(int i=0; i<nparticelle; i++){
                for(int j=i+1; j<nparticelle; j++){
                    double distanza=sqrt((ppvx[nparticelle*i+j]*ppvx[nparticelle*i+j])+(ppvy[nparticelle*i+j]*ppvy[nparticelle*i+j])+(ppvz[nparticelle*i+j]*ppvz[nparticelle*i+j]));
                    if(distanza<lscatola/2.){
                        double r2=(parametripot[1]*parametripot[1])/(distanza*distanza);
                        double r6=r2*r2*r2;
                        double r12=r6*r6;
                        pot1+=r12-r6;
                    }
                }
            }
            double taglio2=(parametripot[1]*parametripot[1])/(lscatola/2.*lscatola/2.);
            double taglio6=taglio2*taglio2*taglio2;
            double taglio12=taglio6*taglio6;
            double pot=parametripot[0]*(pot1-nparticelle*nparticelle/2.*(taglio12-taglio6));
            return pot;
        }
        double logprob(int primomod, int ultimomod){
            /*
            parametri: 
            0) a1
            1) a2
            2) lscatola
            */
            double psi1=0;
            for(int i=primomod;i<ultimomod;i++){
                for(int j=i+1;j<ultimomod;j++){
                    double distanza=sqrt((ppvx[nparticelle*i+j]*ppvx[nparticelle*i+j])+(ppvy[nparticelle*i+j]*ppvy[nparticelle*i+j])+(ppvz[nparticelle*i+j]*ppvz[nparticelle*i+j]));
                    if(distanza<parametrilogprob[2]/2.){
                        psi1+=pow(parametrilogprob[0]/distanza, parametrilogprob[1])+pow(parametrilogprob[0]/(parametrilogprob[2]-distanza), parametrilogprob[1])-2.*pow(parametrilogprob[0]/(parametrilogprob[2]/2.), parametrilogprob[1]);
                    }
                }
            }
            for(int i=0;i<primomod;i++){
                for(int j=primomod;j<ultimomod;j++){
                    double distanza=sqrt((ppvx[nparticelle*i+j]*ppvx[nparticelle*i+j])+(ppvy[nparticelle*i+j]*ppvy[nparticelle*i+j])+(ppvz[nparticelle*i+j]*ppvz[nparticelle*i+j]));
                    if(distanza<parametrilogprob[2]/2.){
                        psi1+=pow(parametrilogprob[0]/distanza, parametrilogprob[1])+pow(parametrilogprob[0]/(parametrilogprob[2]-distanza), parametrilogprob[1])-2.*pow(parametrilogprob[0]/(parametrilogprob[2]/2.), parametrilogprob[1]);
                    }
                }
            }
            for(int i=ultimomod+1;i<nparticelle;i++){
                for(int j=primomod;j<ultimomod;j++){
                    double distanza=sqrt((ppvx[nparticelle*i+j]*ppvx[nparticelle*i+j])+(ppvy[nparticelle*i+j]*ppvy[nparticelle*i+j])+(ppvz[nparticelle*i+j]*ppvz[nparticelle*i+j]));
                    if(distanza<parametrilogprob[2]/2.){
                        psi1+=pow(parametrilogprob[0]/distanza, parametrilogprob[1])+pow(parametrilogprob[0]/(parametrilogprob[2]-distanza), parametrilogprob[1])-2.*pow(parametrilogprob[0]/(parametrilogprob[2]/2.), parametrilogprob[1]);
                    }
                }
            }
            return -psi1;
        }
        int varpv(double variazione[3], int npunto){
            for(int k=0;k<nparticelle;k++){
                ppvx[nparticelle*npunto+k]=(ppvx[nparticelle*npunto+k]+variazione[0])-lscatola*rint((ppvx[nparticelle*npunto+k]+variazione[0])/lscatola);
                ppvx[nparticelle*k+npunto]=(ppvx[nparticelle*k+npunto]-variazione[0])-lscatola*rint((ppvx[nparticelle*k+npunto]-variazione[0])/lscatola);
                ppvy[nparticelle*npunto+k]=(ppvy[nparticelle*npunto+k]+variazione[1])-lscatola*rint((ppvy[nparticelle*npunto+k]+variazione[1])/lscatola);
                ppvy[nparticelle*k+npunto]=(ppvy[nparticelle*k+npunto]-variazione[1])-lscatola*rint((ppvy[nparticelle*k+npunto]-variazione[1])/lscatola);
                ppvz[nparticelle*npunto+k]=(ppvz[nparticelle*npunto+k]+variazione[2])-lscatola*rint((ppvz[nparticelle*npunto+k]+variazione[2])/lscatola);
                ppvz[nparticelle*k+npunto]=(ppvz[nparticelle*k+npunto]-variazione[2])-lscatola*rint((ppvz[nparticelle*k+npunto]-variazione[2])/lscatola);
            }
            return 0;
        }
        Walker(int npart, double lscat, double * pken, int npken, double * ppot, int nppot, double * plogprob, int nplogprob, double * ppvx_gen, double * ppvy_gen, double * ppvz_gen){
            nparticelle=npart; 
            lscatola=lscat;
            peso=1;
            parametriken = new double[npken];
            for(int i=0;i<npken;i++){
                parametriken[i]=pken[i];
            }
            parametripot = new double[nppot];
            for(int i=0;i<nppot;i++){
                parametripot[i]=ppot[i];
            }
            parametrilogprob = new double[nplogprob];
            for(int i=0;i<nplogprob;i++){
                parametrilogprob[i]=plogprob[i];
            }
            ppvx=new double[nparticelle*nparticelle];
            ppvy=new double[nparticelle*nparticelle];
            ppvz=new double[nparticelle*nparticelle];
            for(int i=0;i<nparticelle;i++){
                for(int j=0;j<nparticelle;j++){
                    ppvx[nparticelle*i+j]=ppvx_gen[nparticelle*i+j];
                    ppvy[nparticelle*i+j]=ppvy_gen[nparticelle*i+j];
                    ppvz[nparticelle*i+j]=ppvz_gen[nparticelle*i+j];
                }    
            }
        }
        ~Walker(){
            delete [] ppvx;
            delete [] ppvy;
            delete [] ppvz;
            delete [] parametriken;
            delete [] parametrilogprob;
            delete [] parametripot;
        }
    private:     
        int nparticelle;
        double lscatola=0;
        double * parametriken;
        double * parametripot;
        double * parametrilogprob;
        double * ppvx;
        double * ppvy;
        double * ppvz;
};

int mrt(Walker * walker, double tau, int primomod, int ultimomod, double * (* boxmuller)(double, double, double)){
    int ndamod=ultimomod-primomod+1;
    double variazioni[ndamod][3];
    for(int i=0;i<ndamod;i++){
        double tempA=(double)rand()/(double)RAND_MAX;
        double tempB=(double)rand()/(double)RAND_MAX;
        double tempC=(double)rand()/(double)RAND_MAX;
        double tempD=(double)rand()/(double)RAND_MAX;
        double * tempgaussA=boxmuller(tempA, tempB, sqrt(tau));
        double * tempgaussB=boxmuller(tempC, tempD, sqrt(tau));
        variazioni[i][0]=tempgaussA[0];
        variazioni[i][1]=tempgaussA[1];
        variazioni[i][2]=tempgaussB[0];
        delete [] tempgaussA;
        delete [] tempgaussB;
    }    
    double probprec=(*walker).logprob(primomod, ultimomod);
    for(int i=0;i<ndamod;i++){
        (*walker).varpv(&variazioni[i][0], i);
    }
    double probsucc=(*walker).logprob(primomod, ultimomod);
    double q=exp(-probprec+probsucc);
    if(q>1){
        return 1;
    }
    if((double)rand()/(double)RAND_MAX<q){
        return 1;
    }
    for(int i=0;i<ndamod;i++){
        double menovar[3]={-variazioni[i][0], -variazioni[i][1], -variazioni[i][2]};
        (*walker).varpv(menovar, i);
    }
    return 0;
}
int dmc(Walker walker, double tau, double et, int primomod, int ultimomod, double * (* boxmuller)(double, double, double)){
    int ndamod=ultimomod-primomod+1;
    double variazioni[ndamod][3];
    for(int i=0;i<ndamod;i++){
        double tempA=(double)rand()/(double)RAND_MAX;
        double tempB=(double)rand()/(double)RAND_MAX;
        double tempC=(double)rand()/(double)RAND_MAX;
        double tempD=(double)rand()/(double)RAND_MAX;
        double * tempgaussA=boxmuller(tempA, tempB, sqrt(tau));
        double * tempgaussB=boxmuller(tempC, tempD, sqrt(tau));
        variazioni[i][0]=tempgaussA[0];
        variazioni[i][1]=tempgaussA[1];
        variazioni[i][2]=tempgaussB[0];
        delete [] tempgaussA;
        delete [] tempgaussB;
    }
    double potprec=walker.pot();
    for(int i=0;i<ndamod;i++){
        walker.varpv(&variazioni[i][0], i);
    }
    double potsucc=walker.pot();
    double newpeso=exp(-tau/2.*(potprec+potsucc-2*et));
    walker.peso=walker.peso*0.8+newpeso*0.2;
    return 0;
}

double * mediacumulata(double * vettore, int nelementi, int elementoiniziale){
    double * cumulata;
    cumulata=new double[nelementi-elementoiniziale];
    cumulata[0]=vettore[elementoiniziale];
    for(int i=1;i<nelementi-elementoiniziale;i++){
        cumulata[i]+=cumulata[i-1]+vettore[elementoiniziale+i];
        cumulata[i-1]/=i;
    }
    cumulata[nelementi-elementoiniziale-1]/=nelementi-elementoiniziale;
    return cumulata;
}

double * boxmuller(double x, double y, double sigma){
    double * out;
    out = new double[2];
    out[0]=sigma*sqrt(-2*log(1-y))*cos(2*M_PI*x);
    out[1]=sigma*sqrt(-2*log(1-y))*sin(2*M_PI*x);
    return out;
}

int te_varmc(Walker * walker, int nwork, int * indicafine, int npuntiplot, int nmod, int nparticelle, int * naccettati, double tau, int freqcampionamento, double * potdit, double * kendit, double * stimatore1dit, double * stimatore2dit){
    int contatoreplot=0;
    int npassi=npuntiplot*freqcampionamento;
    for(int i=0; i<npassi; i++){
        int primomod=floor((double)(rand()-1)/(double)RAND_MAX*(nparticelle-nmod));
        int ultimomod=primomod+nmod-1;
        naccettati[nwork]+=mrt(walker, tau, primomod, ultimomod, boxmuller);
        if(i%freqcampionamento==0){
            potdit[npuntiplot*nwork+contatoreplot]=walker->pot();
            double * tempken=walker->ken();
            kendit[npuntiplot*nwork+contatoreplot]=tempken[0];
            stimatore1dit[npuntiplot*nwork+contatoreplot]=tempken[1];
            stimatore2dit[npuntiplot*nwork+contatoreplot]=tempken[2];
            delete [] tempken;
            contatoreplot++;
        }
    }
    *indicafine=1;
    return 0;
}

FILE * valori;

int main(){
    srand(69);
    // parametri che definiscono il tipo di reticolo usato per inizializzare le posizioni
    // int particelleperce=4; // particelle per cella elementare
    // double posizioniince[4][3]={{0, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}}; // posizione della particella nella cella elementare
    int particelleperce=1; // particelle per cella elementare
    double posizioniince[1][3]={{0,0,0}}; // posizione della particella nella cella elementare

    const int nparticelle=64;
    const int npassivar=10000;
    const int freqcampionamento=10;
    const int npuntiplot = npassivar/freqcampionamento+0.5;
    const double densita=0.02;
    const int npassidmc=1000;
    const int nblocchidmc=80;
    const int nwalkers=80;
    const double tau=0.0001;
    const double et=1;
    int nce=(nparticelle/particelleperce);
    double volce=(double)nparticelle/(densita*(double)nce);
    double latoce=pow(volce, 1./3.);
    double latoscatola=latoce*pow(nce, 1./3.);    
    double kendit[nblocchidmc]={0};
    double potdit[nblocchidmc]={0};
    double stimatore1dit[nblocchidmc]={0};
    double stimatore2dit[nblocchidmc]={0};
    double energia[nblocchidmc]={0};
    double kendit_varmc[nwalkers][npuntiplot];
    double potdit_varmc[nwalkers][npuntiplot];
    double stimatore1dit_varmc[nwalkers][npuntiplot];
    double stimatore2dit_varmc[nwalkers][npuntiplot];
    double energiadit_varmc[nwalkers][npuntiplot];
    int naccettati_varmc[nwalkers];
    double * kencum;
    double * potcum;
    double * stimatore1cum;
    double * stimatore2cum;
    double * energiacum;
    double eps=10.22;
    double sigma=2.556;
    double parametripot[2]={eps, sigma};
    double a1=2.5;
    double a2=5;
    double parametriken[2]={a1, a2};
    double parametriprobabilita[4]={a1, a2, latoscatola, (double)nparticelle};
    int nmod=32; // quanti punti modifico ad ogni iterazione del ciclo

    GnuplotDriver kengp;
    GnuplotDriver potgp;
    GnuplotDriver kenditgp;
    GnuplotDriver potditgp;
    GnuplotDriver stimatore1gp;
    GnuplotDriver stimatore2gp;
    GnuplotDriver controllogp;
    GnuplotDriver energiagp;
    GnuplotDriver energiaditgp;
    valori=fopen("dmc/valori.dat", "w");

    // inizializzo le posizioni
    double posx[nparticelle]={0};
    double posy[nparticelle]={0};
    double posz[nparticelle]={0};
    int contatoreparticella=0;
    for(int i=0; i<rint(latoscatola/latoce); i++){
        for(int j=0; j<rint(latoscatola/latoce); j++){
            for(int k=0; k<rint(latoscatola/latoce); k++){
                for(int h=0; h<particelleperce; h++){
                    posx[contatoreparticella]=i*latoce+latoce*posizioniince[h][0];
                    posy[contatoreparticella]=j*latoce+latoce*posizioniince[h][1];
                    posz[contatoreparticella]=k*latoce+latoce*posizioniince[h][2];
                    contatoreparticella++;
                }
            }
        }
    }

    // plot delle posizioni iniziali
    controllogp.conf("fPath", "dmc/posiniziale");
    controllogp.conf("ls", "points");
    controllogp.plot(posx, posy, posz, nparticelle);

    // calcolo ppv
    double * ppvx;
    double * ppvy;
    double * ppvz;
    ppvx = new double[nparticelle*nparticelle];    
    ppvy = new double[nparticelle*nparticelle];
    ppvz = new double[nparticelle*nparticelle];
    for(int i=0; i<nparticelle; i++){
        for(int j=0; j<nparticelle; j++){
            double tempx=posx[j]-posx[i];
            ppvx[i*nparticelle+j]=tempx-latoscatola*rint(tempx/latoscatola);
            double tempy=posy[j]-posy[i];
            ppvy[i*nparticelle+j]=tempy-latoscatola*rint(tempy/latoscatola);
            double tempz=posz[j]-posz[i];
            ppvz[i*nparticelle+j]=tempz-latoscatola*rint(tempz/latoscatola);
        }
    }

    Walker * walkers[nwalkers];
    for(int i=0;i<nwalkers;i++){
        walkers[i]=new Walker(nparticelle, latoscatola, &parametriken[0], 2, &parametripot[0], 2, &parametriprobabilita[0], 4, &ppvx[0], &ppvy[0], &ppvz[0]);
    }
    delete [] ppvx;
    delete [] ppvy;
    delete [] ppvz;

    // evoluzione temporale
        // variational montecarlo
    int nthreads=std::thread::hardware_concurrency()-1;
    std::thread threads[nthreads];
    int nwork=0;
    int finito[nthreads];
    for(int i=0;i<nthreads;i++){
        finito[i]=1;
    }
    while(nwork<nwalkers){
        for(int i=0;i<nthreads;i++){
            if(finito[i]){
                finito[i]=0;
                std::cout<<"\rInizio ad evolvere il walker " << nwork << std::flush;
                threads[i]=std::thread(te_varmc, walkers[nwork], nwork, &finito[i], npuntiplot, nmod, nparticelle, naccettati_varmc, tau, freqcampionamento, &potdit_varmc[0][0], &kendit_varmc[0][0], &stimatore1dit_varmc[0][0], &stimatore2dit_varmc[0][0]);
                threads[i].detach();
                nwork++;
            }
        }
    }
    int completato=0;
    while(!completato){
        completato=1;
        for(int i=0;i<nthreads;i++){
            if(finito[i]==0){
                completato=0;
                break;
            }
        }
    }
    for(int i=0;i<nwalkers;i++){
        std::cout<<"varmc - il walker " << i <<"ha accettazione: " << naccettati_varmc[i] << std::endl;
    }

    potditgp.conf("fPath", "potdmc");
    potditgp.conf("ls", "line");
    potditgp.conf("lim", "[0:1000][-5000:0]");
    for(int i=0;i<nwalkers;i++){
        potditgp.plot(&kendit_varmc[i][0], npuntiplot, 1);
    }

    
    // faccio le medie cumulate

    // grafici

    // calcolo le varianze
    
    // output su file
    
    return 0;
}