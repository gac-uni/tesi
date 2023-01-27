#include "gpd.h"
#include "variables.h"

class Walker{
    public:      
        double * px;
        double * py;
        double * pz;  
        double * ppvx;
        double * ppvy;
        double * ppvz;  
        double * ppvxtemp;
        double * ppvytemp;
        double * ppvztemp;
        double * ptppvx;
        double * ptppvy;
        double * ptppvz;
        double * ptppvxtemp;
        double * ptppvytemp;
        double * ptppvztemp;
        int nparticelle;
        double lscatola;
        double * parametripot;
        double * parametrilogprob;
        double hqsdm = 0.09076;
        double peso=0;
        double interpolapeso=0.8;
        double eT=1;
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
                    double dist=sqrt(ptppvx[nparticelle*j+i]*ptppvx[nparticelle*j+i]+ptppvy[nparticelle*j+i]*ptppvy[nparticelle*j+i]+ptppvz[nparticelle*j+i]*ptppvz[nparticelle*j+i]);
                    if((dist<=lscatola/2.)){
                        lapi+=-1./2.*(dsecondaUT(dist, parametrilogprob)+2.*dprimaUT(dist, parametrilogprob)/dist);
                        gradxi+=1./2.*dprimaUT(dist, parametrilogprob)*ptppvx[nparticelle*j+i]/dist;
                        gradyi+=1./2.*dprimaUT(dist, parametrilogprob)*ptppvy[nparticelle*j+i]/dist;
                        gradzi+=1./2.*dprimaUT(dist, parametrilogprob)*ptppvz[nparticelle*j+i]/dist;
                    }
                }
                for(int j=i+1;j<nparticelle;j++){
                    double dist=sqrt(ptppvx[nparticelle*i+j]*ptppvx[nparticelle*i+j]+ptppvy[nparticelle*i+j]*ptppvy[nparticelle*i+j]+ptppvz[nparticelle*i+j]*ptppvz[nparticelle*i+j]);
                    if((dist<=lscatola/2.)){
                        lapi+=-1./2.*(dsecondaUT(dist, parametrilogprob)+2.*dprimaUT(dist, parametrilogprob)/dist);
                        gradxi+=-1./2.*dprimaUT(dist, parametrilogprob)*ptppvx[nparticelle*i+j]/dist;
                        gradyi+=-1./2.*dprimaUT(dist, parametrilogprob)*ptppvy[nparticelle*i+j]/dist;
                        gradzi+=-1./2.*dprimaUT(dist, parametrilogprob)*ptppvz[nparticelle*i+j]/dist;
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
                    double dist=sqrt(ptppvx[nparticelle*i+j]*ptppvx[nparticelle*i+j]+ptppvy[nparticelle*i+j]*ptppvy[nparticelle*i+j]+ptppvz[nparticelle*i+j]*ptppvz[nparticelle*i+j]);
                    if(dist<=lscatola/2.){
                        pot1+=Vlj(dist, parametripot)-Vlj(lscatola/2., parametripot);
                    }
                }
            }
            return pot1;
        }
        double potdatemp(){
            double pot1=0;
            for(int i=0; i<nparticelle; i++){
                for(int j=i+1; j<nparticelle; j++){
                    double dist=sqrt(ptppvxtemp[nparticelle*i+j]*ptppvxtemp[nparticelle*i+j]+ptppvytemp[nparticelle*i+j]*ptppvytemp[nparticelle*i+j]+ptppvztemp[nparticelle*i+j]*ptppvztemp[nparticelle*i+j]);
                    if(dist<=lscatola/2.){
                        pot1+=Vlj(dist, parametripot) -Vlj(lscatola/2., parametripot);
                    }
                }
            }
            return pot1;
        }
        double logprob(){
            double psi1=0;
            for(int i=0;i<nparticelle;i++){
                for(int j=i+1;j<nparticelle;j++){
                    double dist=sqrt(ptppvx[nparticelle*i+j]*ptppvx[nparticelle*i+j]+ptppvy[nparticelle*i+j]*ptppvy[nparticelle*i+j]+ptppvz[nparticelle*i+j]*ptppvz[nparticelle*i+j]);
                    if(dist<=lscatola/2.){
                        psi1+=UT(dist, parametrilogprob);
                    }
                }
            }
            return -psi1;
        }
        int applicavar(double variazione[3], int npunto){
            px[npunto]+=variazione[0];
            py[npunto]+=variazione[1];
            pz[npunto]+=variazione[2];
            return 0;
        }
        int plotposizioni(std::string nomefile){
            GnuplotDriver posizionigp;
            posizionigp.ls("points");
            posizionigp.limx(-20,20);
            posizionigp.limy(-20,20);
            posizionigp.limz(-20,20);
            posizionigp.fpath(nomefile);
            posizionigp.fext("jpg");
            posizionigp.noprint();
            posizionigp.plot(px, py, pz, nparticelle);
            return 0;
        }
        double davg(){
            double davg=0;
            int count=0;
            for(int i=0;i<nparticelle;i++){
                for(int j=i+1;j<nparticelle;j++){
                    double dist=sqrt(ptppvx[nparticelle*i+j]*ptppvx[nparticelle*i+j]+ptppvy[nparticelle*i+j]*ptppvy[nparticelle*i+j]+ptppvz[nparticelle*i+j]*ptppvz[nparticelle*i+j]);
                    if(dist<=lscatola/2.){
                        davg+=dist;
                        count++;
                    }
                }
            }
            return davg/count;
        }
        double dmin(){
            double dmin=10000;
            int count=0;
            for(int i=0;i<nparticelle;i++){
                for(int j=i+1;j<nparticelle;j++){
                    double dist=sqrt(ptppvx[nparticelle*i+j]*ptppvx[nparticelle*i+j]+ptppvy[nparticelle*i+j]*ptppvy[nparticelle*i+j]+ptppvz[nparticelle*i+j]*ptppvz[nparticelle*i+j]);
                    if(dist<=dmin){
                        dmin=dist;
                    }
                }
            }
            return dmin;
        }
        void computepv(double * locationx, double * locationy, double * locationz){
            for(int i=0;i<nparticelle;i++){
                for(int j=i+1;j<nparticelle;j++){
                    double tempx=px[j]-px[i];
                    double tempy=py[j]-py[i];
                    double tempz=pz[j]-pz[i];
                    locationx[nparticelle*i+j]=tempx-lscatola*rint(tempx/lscatola);
                    locationy[nparticelle*i+j]=tempy-lscatola*rint(tempy/lscatola);
                    locationz[nparticelle*i+j]=tempz-lscatola*rint(tempz/lscatola);
                }
            }
        }
        void computepeso(int accettato, double tau, double * potenziale){
            double vR = pot();
            *potenziale=vR;
            double vRP=0;
            if(accettato){ // allora il passo è stato fatto dallo stato in ptppvxtemp (R') allo stato in ptppvx (R)
                vRP=potdatemp();
            }else{
                vRP=vR;
            }
            double pesoattuale=exp(-tau/2.*(vR+vRP-2*eT));
            peso=interpolapeso*peso+(1.-interpolapeso)*pesoattuale;
        }
        Walker(int npart, double lscat, double * ppot, int nppot, double * plogprob, int nplogprob, double * px_gen, double * py_gen, double * pz_gen){
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
            px=new double[nparticelle];
            py=new double[nparticelle];
            pz=new double[nparticelle];
            for(int j=0;j<nparticelle;j++){
                px[j]=px_gen[j];
                py[j]=py_gen[j];
                pz[j]=pz_gen[j];
            }   
            // inizializzo le matrici dei primi vicini
            ppvx = new double[nparticelle*nparticelle];
            ppvy = new double[nparticelle*nparticelle];
            ppvz = new double[nparticelle*nparticelle];
            ppvxtemp = new double[nparticelle*nparticelle];
            ppvytemp = new double[nparticelle*nparticelle];
            ppvztemp = new double[nparticelle*nparticelle];
            for(int i=0;i<nparticelle;i++){
                for(int j=0;j<nparticelle;j++){
                    ppvx[nparticelle*i+j] = 0;
                    ppvy[nparticelle*i+j] = 0; 
                    ppvz[nparticelle*i+j] = 0; 
                }
            }
            computepv(ppvx, ppvy, ppvz);  
            ptppvx = ppvx;
            ptppvy = ppvy;
            ptppvz = ppvz;
            ptppvxtemp = ppvxtemp;
            ptppvytemp = ppvytemp;
            ptppvztemp = ppvztemp;            
        }        
        ~Walker(){
            delete [] px;
            delete [] py;
            delete [] pz;
            delete [] ppvx;
            delete [] ppvy;
            delete [] ppvz;
            delete [] ppvxtemp;
            delete [] ppvytemp;
            delete [] ppvztemp; 
            delete [] parametrilogprob;
            delete [] parametripot;
        }  
};

class Observables{
    public:
        double * vPot;
        double * vKen;
        double * vStim1;
        double * vStim2;
        double * vEn;
        double * vAcc;
        double * vPesi;
        int nelem;
        Observables(int nelementi){
            nelem = nelementi;
            vPot = new double[nelem];
            vKen = new double[nelem];
            vStim1 = new double[nelem];
            vStim2 = new double[nelem];
            vEn = new double[nelem];
            vAcc = new double[nelem];
            vPesi = new double[nelem];
        }
        ~Observables(){
            delete [] vPot;
            delete [] vKen;
            delete [] vStim1;
            delete [] vStim2;
            delete [] vEn;
            delete [] vAcc;
            delete [] vPesi;
        }
};

double * boxmuller(double x, double y, double sigma){
    double * out;
    out = new double[2];
    out[0]=sigma*sqrt(-2*log(1-y))*cos(2*M_PI*x);
    out[1]=sigma*sqrt(-2*log(1-y))*sin(2*M_PI*x);
    return out;
}

int mrt(Walker * walker, double tau, int count){
    int nparticelle=walker->nparticelle;
    double variazioni[nparticelle][3];
    double probprec=walker->logprob();
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
        double var[3]={variazioni[i][0],variazioni[i][1],variazioni[i][2]};
        walker->applicavar(var, i);
    }
    // sposto i vecchi ppv nei temp
    double * tempx = walker->ptppvxtemp;
    double * tempy = walker->ptppvytemp;
    double * tempz = walker->ptppvztemp;
    walker->ptppvxtemp = walker->ptppvx;
    walker->ptppvytemp = walker->ptppvy;
    walker->ptppvztemp = walker->ptppvz;
    walker->ptppvx=tempx;
    walker->ptppvy=tempy;
    walker->ptppvz=tempz;
    
    // calcolo la configurazione nuova in ppvx, y e z
    walker->computepv(walker->ptppvx,walker->ptppvy,walker->ptppvz);
    double probsucc=walker->logprob();
    double q=exp(-probprec+probsucc); 
    if(q>1){
        return 1;
    }
    if((double)rand()/(double)RAND_MAX<q){
        return 1;
    }
    for(int i=0;i<nparticelle;i++){
        double menovar[3]={-variazioni[i][0],-variazioni[i][1],-variazioni[i][2]};
        walker->applicavar(menovar, i);
    }
    // sposto allo stato originale se rifiuto il passo
    tempx = walker->ptppvxtemp;
    tempy = walker->ptppvytemp;
    tempz = walker->ptppvztemp;
    walker->ptppvxtemp = walker->ptppvx;
    walker->ptppvytemp = walker->ptppvy;
    walker->ptppvztemp = walker->ptppvz;
    walker->ptppvx=tempx;
    walker->ptppvy=tempy;
    walker->ptppvz=tempz;
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

double media(double * v, int nelementi, int elementoiniziale){
    double out = 0;
    for(int i=elementoiniziale;i<nelementi;i++){
        out+=v[i];
    }
    return out/(nelementi-elementoiniziale);
}

void mrt_block(Walker * walker, Observables * obs, int * indicafine){
    int contatoreplot=0;
    double naccettati=0;
    for(int i=0; i<npassiblock; i++){
        int accettato;
        accettato=mrt(walker, tau, i);
        naccettati+=accettato;
        if(i%freqcampionamento==0){
            walker->computepeso(accettato, tau, &(obs->vPot[contatoreplot]));
            obs->vPesi[contatoreplot]=walker->peso;
            double * tempken=walker->ken();
            obs->vKen[contatoreplot]=tempken[0];
            obs->vStim1[contatoreplot]=tempken[1];
            obs->vStim2[contatoreplot]=tempken[2];
            obs->vEn[contatoreplot]=obs->vKen[contatoreplot]+obs->vPot[contatoreplot];
            delete [] tempken;
            contatoreplot++;
            obs->vAcc[contatoreplot]=(double)naccettati/(double)i;
        }
    }
    *indicafine=1;
}

int main(){
    srand(69);

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

    // definizione di walkers e osservabili
    Walker * walkers[nwalkers];
    for(int i=0;i<nwalkers;i++){
        walkers[i] = new Walker(nparticelle, latoscatola, &parametripot[0], 2, &parametrilogprob[0], 2, posx, posy, posz);
    }
    Observables * obs[nwalkers];
    for(int i=0;i<nwalkers;i++){
        obs[i] = new Observables(npassiplot);
    }

    // ciclo sui blocchi di dmc
    std::thread t[n_threads];
    int finito[n_threads];        
    for(int i=0;i<nblocks;i++){ // ciclo sui blocchi di dmc, tra un blocco e l'altro aggiorno eT per mantenere stabile il "numero di walkers"
        // evolvo un blocco di walkers in parallelo
        for(int j=0;j<n_threads;j++){
            finito[j]=1;
        }
        int walkdafare=0;
        while(walkdafare<nwalkers){ // assegno i walkers ai thread finchè ho walkers da assegnare
            for(int j=0;j<n_threads;j++){
                if(finito[j]){
                    finito[j]=0;
                    t[j]=std::thread(mrt_block, walkers[walkdafare], obs[walkdafare], &finito[j]);
                    t[j].detach();
                    walkdafare++;
                    std::cout<<"working on walker: " << walkdafare-1 << std::endl;
                }
                if(walkdafare==nwalkers){
                    break;
                }
            }
        }
        bool completato=false; // devo aspettare che tutti i thread finiscano
        while(!completato){
            completato=true;
            for(int j=0;j<n_threads;j++){
                if(!finito[j]){
                    completato=false;
                    break;
                }
            }
        }
    }

    // plotto il potenziale, le energie cinetiche ed i pesi per i vari walkers
    GnuplotDriver kengp;
    GnuplotDriver stim1gp;
    GnuplotDriver stim2gp;
    GnuplotDriver potgp;
    GnuplotDriver pesigp;

    kengp.fpath("dmc_out/ken");
    kengp.ls("line");
    stim1gp.fpath("dmc_out/stim1");
    stim1gp.ls("line");
    stim2gp.fpath("dmc_out/stim2");
    stim2gp.ls("line");
    potgp.fpath("dmc_out/pot");
    potgp.ls("line");
    pesigp.fpath("dmc_out/pesi");
    pesigp.ls("line");

    for(int i=0;i<nwalkers;i++){
        kengp.plot(obs[i]->vKen, npassiplot, 1);
        stim1gp.plot(obs[i]->vStim1, npassiplot, 1);
        stim2gp.plot(obs[i]->vStim2, npassiplot, 1);
        potgp.plot(obs[i]->vPot, npassiplot, 1);
        pesigp.plot(obs[i]->vPesi, npassiplot, 1);
    }

    // calcolo i valori numerici di medie e medie dei quadrati pesate con i pesi di dmc
    double ken[npassiplot];
    double kenquad[npassiplot];
    double stim1[npassiplot];
    double stim1quad[npassiplot];
    double stim2[npassiplot];
    double stim2quad[npassiplot];
    double pot[npassiplot];
    double potquad[npassiplot];
    for(int i=0;i<npassiplot;i++){
        double kenloc=0;
        double kenquadloc=0;
        double stim1loc=0;
        double stim1quadloc=0;
        double stim2loc=0;
        double stim2quadloc=0;
        double potloc=0;
        double potquadloc=0;
        double pesiloc=0;
        for(int j=0;j<nwalkers;j++){
            kenloc+=obs[j]->vKen[i]*obs[j]->vPesi[i];
            kenquadloc+=obs[j]->vKen[i]*obs[j]->vKen[i]*obs[j]->vPesi[i];
            stim1loc+=obs[j]->vStim1[i]*obs[j]->vPesi[i];
            stim1quadloc+=obs[j]->vStim1[i]*obs[j]->vStim1[i]*obs[j]->vPesi[i];
            stim2loc+=obs[j]->vStim2[i]*obs[j]->vPesi[i];
            stim2quadloc+=obs[j]->vStim2[i]*obs[j]->vStim2[i]*obs[j]->vPesi[i];
            potloc+=obs[j]->vPot[i]*obs[j]->vPesi[i];
            potquadloc+=obs[j]->vPot[i]*obs[j]->vPot[i]*obs[j]->vPesi[i];
            pesiloc+=obs[j]->vPesi[i];            
        }
        ken[i] = kenloc/pesiloc;
        kenquad[i] = kenquadloc/pesiloc;
        stim1[i] = stim1loc/pesiloc;
        stim1quad[i] = stim1quadloc/pesiloc;
        stim2[i] = stim2loc/pesiloc;
        stim2quad[i] = stim2quadloc/pesiloc;
        pot[i] = potloc/pesiloc;
        potquad[i] = potquadloc/pesiloc;
    }

    // plotto i valori medi tra le varie run di dmc
    GnuplotDriver kenpesgp;
    GnuplotDriver stim1pesgp;
    GnuplotDriver stim2pesgp;
    GnuplotDriver potpesgp;

    kenpesgp.fpath("dmc_out/kenpes");
    kenpesgp.ls("line");
    kenpesgp.plot(ken, npassiplot, 1);
    stim1pesgp.fpath("dmc_out/stim1pes");
    stim1pesgp.ls("line");
    stim1pesgp.plot(stim1, npassiplot, 1);
    stim2pesgp.fpath("dmc_out/stim2pes");
    stim2pesgp.ls("line");
    stim2pesgp.plot(stim2, npassiplot, 1);
    potpesgp.fpath("dmc_out/potpes");
    potpesgp.ls("line");
    potpesgp.plot(pot, npassiplot, 1);

    // calcolo i valori medi totali e le varianze e le stampo su file
    double kenmed=media(ken, npassiplot, 100);
    double kenquadmed=media(kenquad, npassiplot, 100);
    double stim1med=media(stim1, npassiplot, 100);
    double stim1quadmed=media(stim1quad, npassiplot, 100);
    double stim2med=media(stim2, npassiplot, 100);
    double stim2quadmed=media(stim2quad, npassiplot, 100);
    double potmed=media(pot, npassiplot, 100);
    double potquadmed=media(potquad, npassiplot, 100);

    fprintf(valori, "ken medio: %lf con varianza: %lf\n", kenmed, kenquadmed-kenmed*kenmed);
    fprintf(valori, "stim1 medio: %lf con varianza: %lf\n", stim1med, stim1quadmed-stim1med*stim1med);
    fprintf(valori, "stim2 medio: %lf con varianza: %lf\n", stim2med, stim2quadmed-stim2med*stim2med);
    fprintf(valori, "pot medio: %lf con varianza: %lf\n", potmed, potquadmed-potmed*potmed);

    return 0;
}