#include <iostream>
#include <cmath>
#include <cstdio>
#include <cassert>
#include "gpd.h"

class Walker{
    public:      
        double * ppvx;
        double * ppvy;
        double * ppvz;
        int nparticelle;
        double lscatola;
        double * parametripot;
        double * parametrilogprob;
        double hqsdm = 0.09076;
        double dzeroU(double r, double * parametri){
            /*
            parametri: 
            0) a1
            1) a2
            */
            return pow(parametri[0]/r, parametri[1]);
        }
        double dprimaU(double r, double * parametri){
            /*
            parametri: 
            0) a1
            1) a2
            */
            return -pow(parametri[0], parametri[1])*parametri[1]*pow(r, -parametri[1]-1);
        }
        double dsecondaU(double r, double * parametri){
            /*
            parametri: 
            0) a1
            1) a2
            */
            return pow(parametri[0], parametri[1])*parametri[1]*(parametri[1]+1)*pow(r, -parametri[1]-2);
        }
        double UT(double r, double * parametri){
            /*
            parametri: 
            0) a1
            1) a2
            */
            return dzeroU(r, parametri)+dzeroU(lscatola-r, parametri)-2*dzeroU(lscatola/2., parametri);
        }
        double dprimaUT(double r, double * parametri){
            /*
            parametri: 
            0) a1
            1) a2
            */
            return dprimaU(r, parametri)-dprimaU(lscatola-r, parametri);
        }
        double dsecondaUT(double r, double * parametri){
            /*
            parametri: 
            0) a1
            1) a2
            */
            return dsecondaU(r, parametri)+dsecondaU(lscatola-r, parametri);
        }
        double * ken(){
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
            double laptot=0;
            double gradquadtot=0;
            for(int i=0;i<nparticelle; i++){
                double lapi=0;
                double gradxi=0;
                double gradyi=0;
                double gradzi=0;
                for(int j=0;j<i;j++){                    
                    double dist=sqrt((ppvx[nparticelle*i+j]*ppvx[nparticelle*i+j])+(ppvy[nparticelle*i+j]*ppvy[nparticelle*i+j])+(ppvz[nparticelle*i+j]*ppvz[nparticelle*i+j]));
                    if((dist<=lscatola/2.)){
                        lapi+=+-1./2.*(dsecondaUT(dist, parametrilogprob)+2.*dprimaUT(dist, parametrilogprob)/dist);
                        gradxi+=-1./2.*dprimaUT(dist, parametrilogprob)*ppvx[nparticelle*i+j]/dist;
                        gradyi+=-1./2.*dprimaUT(dist, parametrilogprob)*ppvy[nparticelle*i+j]/dist;
                        gradzi+=-1./2.*dprimaUT(dist, parametrilogprob)*ppvz[nparticelle*i+j]/dist;
                    }
                }
                for(int j=i+1;j<nparticelle;j++){
                    double dist=sqrt((ppvx[nparticelle*i+j]*ppvx[nparticelle*i+j])+(ppvy[nparticelle*i+j]*ppvy[nparticelle*i+j])+(ppvz[nparticelle*i+j]*ppvz[nparticelle*i+j]));
                    if((dist<=lscatola/2.)){
                        lapi+=-1./2.*(dsecondaUT(dist, parametrilogprob)+2.*dprimaUT(dist, parametrilogprob)/dist);
                        gradxi+=-1./2.*dprimaUT(dist, parametrilogprob)*ppvx[nparticelle*i+j]/dist;
                        gradyi+=-1./2.*dprimaUT(dist, parametrilogprob)*ppvy[nparticelle*i+j]/dist;
                        gradzi+=-1./2.*dprimaUT(dist, parametrilogprob)*ppvz[nparticelle*i+j]/dist;
                    }
                }
                laptot+=lapi;
                gradquadtot+=gradxi*gradxi+gradyi*gradyi+gradzi*gradzi;
            }
            out[0]=-hqsdm*(laptot+gradquadtot);
            out[1]=hqsdm*gradquadtot;
            out[2]=-1./2.*hqsdm*laptot;
            return out;
        }
        double Vlj(double r, double * parametri){
            double alla3=(parametri[1]/r)*(parametri[1]/r)*(parametri[1]/r);
            double alla6=alla3*alla3;
            double alla12=alla6*alla6;
            return 4*parametri[0]*(alla12-alla6);
        }
        double pot(){
            double pot1=0;
            for(int i=0; i<nparticelle; i++){
                for(int j=i+1; j<nparticelle; j++){                       
                    double dist=sqrt((ppvx[nparticelle*i+j]*ppvx[nparticelle*i+j])+(ppvy[nparticelle*i+j]*ppvy[nparticelle*i+j])+(ppvz[nparticelle*i+j]*ppvz[nparticelle*i+j]));
                    if(dist<=lscatola/2.){
                        pot1+=Vlj(dist, parametripot)-Vlj(lscatola/2., parametripot);
                    }
                }
            }
            return pot1;
        }
        double logprob(){
            double psi1=0;
            for(int i=0;i<nparticelle;i++){
                for(int j=i+1;j<nparticelle;j++){
                    double distanza=sqrt((ppvx[nparticelle*i+j]*ppvx[nparticelle*i+j])+(ppvy[nparticelle*i+j]*ppvy[nparticelle*i+j])+(ppvz[nparticelle*i+j]*ppvz[nparticelle*i+j]));
                    if(distanza<=lscatola/2.){
                        psi1+=UT(distanza, parametrilogprob);
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
        int plotposizioni(std::string nomefile){
            GnuplotDriver gp;
            gp.fpath(nomefile);
            gp.fext("jpg");
            gp.ls("points");
            gp.noprint();
            gp.plot(ppvx, ppvy, ppvz, nparticelle);
            return 0;
        }
        Walker(int npart, double lscat, double * ppot, int nppot, double * plogprob, int nplogprob, double * ppvx_gen, double * ppvy_gen, double * ppvz_gen){
            nparticelle=npart; 
            lscatola=lscat;
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
            delete [] parametrilogprob;
            delete [] parametripot;
        }  
};

double * boxmuller(double x, double y, double sigma){
    double * out;
    out = new double[2];
    out[0]=sigma*sqrt(-2*log(1-y))*cos(2*M_PI*x);
    out[1]=sigma*sqrt(-2*log(1-y))*sin(2*M_PI*x);
    return out;
}

int mrt(Walker * walker, double tau){
    int nparticelle = walker->nparticelle;
    double variazioni[nparticelle][3];
    double probprec=(*walker).logprob();
    for(int i=0;i<nparticelle;i++){
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
    for(int i=0;i<nparticelle;i++){
        (*walker).varpv(&variazioni[i][0], i);
    }
    double probsucc=(*walker).logprob();
    double q=exp(-probprec+probsucc); 
    if(q>1){
        return 1;
    }
    if((double)rand()/(double)RAND_MAX<q){
        return 1;
    }
    for(int i=0;i<nparticelle;i++){
        double menovar[3]={-variazioni[i][0], -variazioni[i][1], -variazioni[i][2]};
        (*walker).varpv(menovar, i);
    }
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

FILE * valori;
FILE * out;

int main(){
    srand(69);
    out = fopen("out.dat", "w");
    // parametri che definiscono il tipo di reticolo usato per inizializzare le posizioni
    // int particelleperce=4; // particelle per cella elementare
    // double posizioniince[4][3]={{0, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}}; // posizione della particella nella cella elementare
    int particelleperce=1; // particelle per cella elementare
    double posizioniince[1][3]={{0,0,0}}; // posizione della particella nella cella elementare

    int nparticelle=64;
    double eps=10.22;
    double sigma=2.556;
    double densita=0.0218*sigma*sigma*sigma; // in unità adimensionali \tilda(r)=r/sigma
    int nce=(nparticelle/particelleperce);
    double volce=(double)nparticelle/(densita*(double)nce);
    double latoce=pow(volce, 1./3.);
    double latoscatola=latoce*pow(nce, 1./3.);
    int npassi=300000;
    int freqcampionamento=300;
    int npassiplot = npassi/freqcampionamento+0.5; 
    double tau=0.0008;
    double naccettati=0;
    double kendit[npassiplot]={0};
    double potdit[npassiplot]={0};
    double stimatore1dit[npassiplot]={0};
    double stimatore2dit[npassiplot]={0};
    double energia[npassiplot]={0};
    double * kencum;
    double * potcum;
    double * stimatore1cum;
    double * stimatore2cum;
    double * energiacum;
    double parametripot[2]={1, 1};
    double a1=0.978;
    double a2=5;
    double parametrilogprob[2]={a1, a2};

    GnuplotDriver kengp;
    GnuplotDriver potgp;
    GnuplotDriver kenditgp;
    GnuplotDriver potditgp;
    GnuplotDriver stimatore1gp;
    GnuplotDriver stimatore2gp;
    GnuplotDriver controllogp;
    GnuplotDriver energiagp;
    GnuplotDriver energiaditgp;
    valori=fopen("varmcbello_out/valori.dat", "w");

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
    controllogp.fpath("varmcbello_out/posiniziale");
    controllogp.ls("points");
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

    Walker walker(nparticelle, latoscatola, &parametripot[0], 2, &parametrilogprob[0], 2, &ppvx[0], &ppvy[0], &ppvz[0]);
    delete [] ppvx;
    delete [] ppvy;
    delete [] ppvz;

    // evoluzione temporale
    int contatoreplot=0;
    double accettazione[(int)rint(npassi/freqcampionamento)]={0};
    for(int i=0; i<npassi; i++){
        naccettati+=mrt(&walker, tau);
        if(i%freqcampionamento==0){
            potdit[contatoreplot]=walker.pot();
            double * tempken=walker.ken();
            kendit[contatoreplot]=tempken[0];
            stimatore1dit[contatoreplot]=tempken[1];
            stimatore2dit[contatoreplot]=tempken[2];
            energia[contatoreplot]=kendit[contatoreplot]+potdit[contatoreplot];
            delete [] tempken;
            contatoreplot++;
            accettazione[contatoreplot]=(double)naccettati/(double)i;
            // salvo le posizioni dei punti come plot
            std::string nome = "varmcbello_out/pos/pos";
            nome = nome+std::to_string(contatoreplot);
            walker.plotposizioni(nome);
            std::cout<<"\rcompletamento: " << (double)i/(double)npassi*100. << "%, probabilità di accettazione: " << (double)naccettati/(double)i << std::flush;
        }
    }
    std::cout<<std::endl;

    // faccio le medie cumulate
    kencum=mediacumulata(kendit, npassiplot, 0);
    potcum=mediacumulata(potdit, npassiplot, 0);
    stimatore1cum=mediacumulata(stimatore1dit, npassiplot, 0);
    stimatore2cum=mediacumulata(stimatore2dit, npassiplot, 0);
    energiacum=mediacumulata(energia, npassiplot, 0);

    for(int i=0;i<npassiplot; i++){
        fprintf(out, "%lf %lf\n", kendit[i], potdit[i]);
    }

    // grafici
    // kengp.fpath("varmcbello_out/ken");
    // kengp.ls("line");
    // kengp.plot(kencum, npassiplot, 1);

    // potgp.fpath("varmcbello_out/pot");
    // potgp.ls("line");
    // potgp.plot(potcum, npassiplot, 1);

    kenditgp.fpath("varmcbello_out/kendit");
    kenditgp.ls("line");
    kenditgp.t("energia cinetica");
    kenditgp.y("energia [adimensionale]");
    kenditgp.plot(kendit, npassiplot, 1);

    potditgp.fpath("varmcbello_out/potdit");
    potditgp.ls("line");
    potditgp.t("potenziale");
    potditgp.y("energia [adimensionale]");
    potditgp.plot(potdit, npassiplot, 1);

    stimatore1gp.fpath("varmcbello_out/stimatore1");
    stimatore1gp.ls("line");
    stimatore1gp.plot(stimatore1dit, npassiplot, 1);

    stimatore1gp.fpath("varmcbello_out/stimatore2");
    stimatore1gp.ls("line");
    stimatore2gp.plot(stimatore2dit, npassiplot, 1);

    // energiagp.fpath("varmcbello_out/energia");
    // energiagp.ls("line");
    // energiagp.plot(energiacum, npassiplot, 1)

    energiaditgp.fpath("varmcbello_out/energiadit");
    energiaditgp.ls("line");
    energiaditgp.t("potenziale");
    energiaditgp.y("energia [adimensionale]");
    energiaditgp.plot(energia, npassiplot, 1);

    // calcolo le varianze
    double potquadmedio=0;
    double kenquadmedio=0;
    double stim1quadmedio=0;
    double stim2quadmedio=0;
    double energiaquadmedio=0;
    for(int i=0;i<npassiplot;i++){
        potquadmedio+=potdit[i]*potdit[i]/npassi;
        kenquadmedio+=kendit[i]*kendit[i]/npassi;
        stim1quadmedio+=stimatore1dit[i]*stimatore1dit[i]/npassi;
        stim2quadmedio+=stimatore2dit[i]*stimatore2dit[i]/npassi;
        energiaquadmedio+=energia[i]*energia[i];
    }
    double varken=kenquadmedio-(kencum[npassiplot-1]*kencum[npassiplot-1]);
    double varpot=potquadmedio-(potcum[npassiplot-1]*potcum[npassiplot-1]);
    double varstim1=stim1quadmedio-(stimatore1cum[npassiplot-1]*stimatore1cum[npassiplot-1]);
    double varstim2=stim2quadmedio-(stimatore2cum[npassiplot-1]*stimatore2cum[npassiplot-1]);
    double varenergia=energiaquadmedio-(energiacum[npassiplot-1]*energiacum[npassiplot-1]);

    // output su file
    fprintf(valori, "energia: %lf, con varianza: %lf\n", energiacum[npassiplot-1], varenergia);
    fprintf(valori, "energia cinetica: %lf, con varianza: %lf\n", kencum[npassiplot-1], varken);
    fprintf(valori, "potenziale: %lf, con varianza: %lf\n", potcum[npassiplot-1], varpot);
    fprintf(valori, "primo stimatore: %lf, con varianza: %lf\n", stimatore1cum[npassiplot-1], varstim1);
    fprintf(valori, "secondo stimatore: %lf, con varianza: %lf\n", stimatore2cum[npassiplot-1], varstim2);

    return 0;
}