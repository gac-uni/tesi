#include <iostream>
#include <cassert>
#include "variables_mt.h"

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
        int computegradi(int i, double * gradx, double * grady, double * gradz){
            *gradx=0;
            *grady=0;
            *gradz=0;
            for(int j=0;j<i;j++){
                double dist=sqrt(ptppvx[nparticelle*j+i]*ptppvx[nparticelle*j+i]+ptppvy[nparticelle*j+i]*ptppvy[nparticelle*j+i]+ptppvz[nparticelle*j+i]*ptppvz[nparticelle*j+i]);
                if((dist<=lscatola/2.)){
                    *gradx+= -1./2.*dprimaUT(dist, parametrilogprob)*ptppvx[nparticelle*j+i]/dist;
                    *grady+= -1./2.*dprimaUT(dist, parametrilogprob)*ptppvy[nparticelle*j+i]/dist;
                    *gradz+= -1./2.*dprimaUT(dist, parametrilogprob)*ptppvz[nparticelle*j+i]/dist;
                }
            }
            for(int j=i+1;j<nparticelle;j++){
                double dist=sqrt(ptppvx[nparticelle*i+j]*ptppvx[nparticelle*i+j]+ptppvy[nparticelle*i+j]*ptppvy[nparticelle*i+j]+ptppvz[nparticelle*i+j]*ptppvz[nparticelle*i+j]);
                if((dist<=lscatola/2.)){
                    *gradx+= 1./2.*dprimaUT(dist, parametrilogprob)*ptppvx[nparticelle*i+j]/dist;
                    *grady+= 1./2.*dprimaUT(dist, parametrilogprob)*ptppvy[nparticelle*i+j]/dist;
                    *gradz+= 1./2.*dprimaUT(dist, parametrilogprob)*ptppvz[nparticelle*i+j]/dist;
                }
            }
            return 0;
        }
        double * ken(double * gradlogattx, double * gradlogatty, double * gradlogattz){
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
                    }
                }
                for(int j=i+1;j<nparticelle;j++){
                    double dist=sqrt(ptppvx[nparticelle*i+j]*ptppvx[nparticelle*i+j]+ptppvy[nparticelle*i+j]*ptppvy[nparticelle*i+j]+ptppvz[nparticelle*i+j]*ptppvz[nparticelle*i+j]);
                    if((dist<=lscatola/2.)){
                        lapi+=-1./2.*(dsecondaUT(dist, parametrilogprob)+2.*dprimaUT(dist, parametrilogprob)/dist);
                    }
                }
                // computegradi(i, &gradxi, &gradyi, &gradzi);
                // gradlogprex[i] = gradxi;
                // gradlogprey[i] = gradyi;
                // gradlogprez[i] = gradzi;
                gradxi=gradlogattx[i];
                gradyi=gradlogatty[i];
                gradzi=gradlogattz[i];
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
        // int plotposizioni(std::string nomefile){
        //     GnuplotDriver posizionigp;
        //     posizionigp.ls("points");
        //     posizionigp.limx(-20,20);
        //     posizionigp.limy(-20,20);
        //     posizionigp.limz(-20,20);
        //     posizionigp.fpath(nomefile);
        //     posizionigp.fext("jpg");
        //     posizionigp.plot(px, py, pz, nparticelle);
        //     return 0;
        // }
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
        // double * vPesi;
        int nelem;
        Observables(int nelementi){
            nelem = nelementi;
            vPot = new double[nelem];
            vKen = new double[nelem];
            vStim1 = new double[nelem];
            vStim2 = new double[nelem];
            vEn = new double[nelem];
            vAcc = new double[nelem];
            // vPesi = new double[nelem];
        }
        ~Observables(){
            delete [] vPot;
            delete [] vKen;
            delete [] vStim1;
            delete [] vStim2;
            delete [] vEn;
            delete [] vAcc;
            // delete [] vPesi;
        }
};

double * boxmuller(double x, double y, double sigma){
    double * out;
    out = new double[2];
    out[0]=sigma*sqrt(-2*log(1-y))*cos(2*M_PI*x);
    out[1]=sigma*sqrt(-2*log(1-y))*sin(2*M_PI*x);
    return out;
}

int montecarlo(Walker * walker, double tau, double * gradlnpsiprex, double * gradlnpsiprey, double * gradlnpsiprez, double * gradlnpsisuccx, double * gradlnpsisuccy, double * gradlnpsisuccz){
    int nparticelle=walker->nparticelle;
    double variazioni[nparticelle][3];
    double psiprec=walker->logprob();
    for(int i=0;i<nparticelle;i++){
        double tempA=(double)rand()/(double)RAND_MAX;
        double tempB=(double)rand()/(double)RAND_MAX;
        double tempC=(double)rand()/(double)RAND_MAX;
        double tempD=(double)rand()/(double)RAND_MAX;
        double * tempgaussA=boxmuller(tempA, tempB, sqrt(tau));
        double * tempgaussB=boxmuller(tempC, tempD, sqrt(tau));
        variazioni[i][0]=tempgaussA[0]+tau*gradlnpsiprex[i];
        variazioni[i][1]=tempgaussA[1]+tau*gradlnpsiprey[i];
        variazioni[i][2]=tempgaussB[0]+tau*gradlnpsiprez[i];
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

    // calcolo i gradienti per la drift velocity nella nuova configurazione
    for(int i=0;i<nparticelle;i++){
        walker->computegradi(i, &gradlnpsisuccx[i], &gradlnpsisuccy[i], &gradlnpsisuccz[i]);
    }
    double psisucc=walker->logprob();

    // calcolo i prodotti scalari necessari a calcolare le due G_d (eq 3.49 e 3.52)
    // calcolo tutto considerando le coordinate come le coordinate dei primi vicini centrate nella prima particella perche mi viene più comodo ma potrebbe essere sbagliato
    double exp_g_d_dritta=0;
    double exp_g_d_storta=0;
    for(int i=0;i<nparticelle;i++){
        double posxpre=walker->px[i]-variazioni[i][0];
        double posypre=walker->py[i]-variazioni[i][1];
        double poszpre=walker->pz[i]-variazioni[i][2];
        exp_g_d_dritta += -(walker->px[i]-posxpre-tau*gradlnpsiprex[i])*(walker->px[i]-posxpre-tau*gradlnpsiprex[i]);
        exp_g_d_dritta += -(walker->py[i]-posypre-tau*gradlnpsiprey[i])*(walker->py[i]-posypre-tau*gradlnpsiprey[i]);
        exp_g_d_dritta += -(walker->pz[i]-poszpre-tau*gradlnpsiprez[i])*(walker->pz[i]-poszpre-tau*gradlnpsiprez[i]);
        exp_g_d_storta += -(posxpre-walker->px[i]-tau*gradlnpsisuccx[i])*(posxpre-walker->px[i]-tau*gradlnpsisuccx[i]);
        exp_g_d_storta += -(posypre-walker->py[i]-tau*gradlnpsisuccy[i])*(posypre-walker->py[i]-tau*gradlnpsisuccy[i]);
        exp_g_d_storta += -(poszpre-walker->pz[i]-tau*gradlnpsisuccz[i])*(poszpre-walker->pz[i]-tau*gradlnpsisuccz[i]);
    }
    exp_g_d_dritta/=2.*tau;
    exp_g_d_storta/=2.*tau;

    // posso calcolare la probabilità di accettazione
    double q=exp(-exp_g_d_dritta+exp_g_d_storta-psiprec+psisucc); 
    if(q>1){
        return 1;
    }
    if((double)rand()/(double)RAND_MAX<q){
        return 1;
    }
    // se rifiuto il passo sistemo le posizioni allo stato precedente
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
    // se rifiuto il passo i gradienti successivi sono uguali ai gradienti precedenti
    for(int i=0;i<nparticelle;i++){
        gradlnpsisuccx[i] = gradlnpsiprex[i];
        gradlnpsisuccy[i] = gradlnpsiprey[i];
        gradlnpsisuccz[i] = gradlnpsiprez[i]; 
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

double * media(double * v, int nelementi, int elementoiniziale){
    double * out = new double[2];
    out[0]=0; 
    double medquad=0;
    for(int i=elementoiniziale;i<nelementi;i++){
        out[0]+=v[i];
        medquad +=v[i]*v[i];
    }
    out[0]/=(nelementi-elementoiniziale);
    medquad/=(nelementi-elementoiniziale);
    out[1] = medquad-out[0]*out[0];
    return out;
}

void mrt_block(Walker * walker, Observables * obs, int * indicafine){
    int contatoreplot=0;
    double naccettati=0;
    double * gradlogprex;
    gradlogprex = new double[nparticelle];
    double * gradlogprey;
    gradlogprey = new double[nparticelle];
    double * gradlogprez;
    gradlogprez = new double[nparticelle];
    double * gradlogsuccx;
    gradlogsuccx = new double[nparticelle];
    double * gradlogsuccy;
    gradlogsuccy = new double[nparticelle];
    double * gradlogsuccz;
    gradlogsuccz = new double[nparticelle];
    for(int i=0; i<npassiblock; i++){
        if(i==0){
            for(int j=0;j<nparticelle;j++){
                walker->computegradi(j, &gradlogprex[j], &gradlogprey[j], &gradlogprez[j]);
            }
        }
        int accettato;
        accettato=montecarlo(walker, tau, gradlogprex, gradlogprey, gradlogprez, gradlogsuccx, gradlogsuccy, gradlogsuccz);
        naccettati+=accettato;
        // if(i%freqcampionamento==0){
        double * tempken=walker->ken(gradlogsuccx, gradlogsuccy, gradlogsuccz);
        obs->vKen[contatoreplot]=tempken[0];
        obs->vStim1[contatoreplot]=tempken[1];
        obs->vStim2[contatoreplot]=tempken[2];
        obs->vPot[contatoreplot]=walker->pot();
        obs->vEn[contatoreplot]=obs->vKen[contatoreplot]+obs->vPot[contatoreplot];
        delete [] tempken;
        contatoreplot++;
        obs->vAcc[contatoreplot]=(double)naccettati/(double)i;
        for(int j=0;j<nparticelle;j++){
            gradlogprex[j]=gradlogsuccx[j];
            gradlogprey[j]=gradlogsuccy[j];
            gradlogprez[j]=gradlogsuccz[j];
        }
        // }
    }
    delete [] gradlogprex;
    delete [] gradlogprey;
    delete [] gradlogprez;
    delete [] gradlogsuccx;
    delete [] gradlogsuccy;
    delete [] gradlogsuccz;
    

    *indicafine=1;
}

int main(){
    // srand(69);


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
    std::vector<Walker *> walkers(nwalkers);
    for(int i=0;i<nwalkers;i++){
        walkers[i] = new Walker(nparticelle, latoscatola, &parametripot[0], 2, &parametrilogprob[0], 2, posx, posy, posz);
    }
    std::vector<Observables *> obs(nwalkers);
    for(int i=0;i<nwalkers;i++){
        obs[i] = new Observables(npassiblock);
    }
    double et[npassidmc];
    std::vector<double> pesi(nwalkers);
    for(int i=0;i<nwalkers;i++){
        pesi[i]=1;
    }
    double pwanted = nwalkers;
    double * energiemediebig = new double[npassidmc*npassiblock];
    for(int i=0;i<npassidmc*npassiblock;i++){
        energiemediebig[i]=0;
    }
    double * outpesi = new double[npassidmc];
    for(int i=0;i<npassidmc;i++){
        outpesi[i]=0;
    }
    double * stim2out = new double[npassidmc*npassiblock];
    for(int i=0;i<npassidmc*npassiblock;i++){
        stim2out[i]=0;
    }
    double * kenout = new double[npassidmc*npassiblock];
    for(int i=0;i<npassidmc*npassiblock;i++){
        kenout[i]=0;
    }
    double * potout = new double[npassidmc*npassiblock];
    for(int i=0;i<npassidmc*npassiblock;i++){
        potout[i]=0;
    }
    double * stim1out = new double[npassidmc*npassiblock];
    for(int i=0;i<npassidmc*npassiblock;i++){
        stim1out[i]=0;
    }
    for(int countdmc=0;countdmc<npassidmc;countdmc++){
        // std::cout << "blocco nr: " << countdmc << std::endl << std::endl;
        // ciclo sui blocchi
        std::thread t[n_threads];
        int finito[n_threads];        
        // for(int i=0;i<nblocks;i++){ // ciclo sui blocchi di dmc, tra un blocco e l'altro aggiorno eT per mantenere stabile il "numero di walkers"
        // evolvo un blocco di walkers in parallelo
        for(int j=0;j<n_threads;j++){
            finito[j]=1;
        }
        int walkdafare=0;
        // std::cout << countdmc << " " << nwalkers << std::endl;
        while(walkdafare<nwalkers){ // assegno i walkers ai thread finchè ho walkers da assegnare
            for(int j=0;j<n_threads;j++){
                if(finito[j]){
                    finito[j]=0;
                    t[j]=std::thread(mrt_block, walkers[walkdafare], obs[walkdafare], &finito[j]);
                    t[j].detach();
                    walkdafare++;
                    // std::cout<<"working on walker: " << walkdafare-1 << std::endl;
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
        // calcolo le energie medie
        for(int i=0;i<npassiblock;i++){
            for(int j=0;j<nwalkers;j++){
                energiemediebig[countdmc*npassiblock+i] += obs[j]->vKen[i] + obs[j]->vPot[i]; 
            }
            energiemediebig[countdmc*npassiblock+i]/=nwalkers;
        }
        // calcolo stimatori medi
        for(int i=0;i<npassiblock;i++){
            for(int j=0;j<nwalkers;j++){
                kenout[countdmc*npassiblock+i] += obs[j]->vKen[i]; 
            }
            kenout[countdmc*npassiblock+i]/=nwalkers;
        }
        for(int i=0;i<npassiblock;i++){
            for(int j=0;j<nwalkers;j++){
                potout[countdmc*npassiblock+i] += obs[j]->vPot[i]; 
            }
            potout[countdmc*npassiblock+i]/=nwalkers;
        }
        for(int i=0;i<npassiblock;i++){
            for(int j=0;j<nwalkers;j++){
                stim1out[countdmc*npassiblock+i] += obs[j]->vStim1[i]; 
            }
            stim1out[countdmc*npassiblock+i]/=nwalkers;
        }
        for(int i=0;i<npassiblock;i++){
            for(int j=0;j<nwalkers;j++){
                stim2out[countdmc*npassiblock+i] += obs[j]->vStim2[i]; 
            }
            stim2out[countdmc*npassiblock+i]/=nwalkers;
        }

        //calcolo i de
        double detemp[nwalkers]={0};
        for(int i=0;i<nwalkers;i++){
            double kenprec;
            double kensucc;
            double potprec;
            double potsucc;
            kensucc=obs[i]->vKen[0];
            potsucc=obs[i]->vPot[0];
            for(int j=0;j<npassiblock-1;j++){
                kenprec=kensucc;
                potprec=potsucc;
                kensucc=obs[i]->vKen[j+1];
                potsucc=obs[i]->vPot[j+1];
                detemp[i]+=kenprec+potprec+kensucc+potsucc;
            }
        }
        //calcolo et 
        double sommaperet=0;
        for(int i=0;i<nwalkers;i++){
            sommaperet+=pesi[i]*exp(-tau/2.*detemp[i]);
        }
        et[countdmc]=1/(tau*npassidmc)*(log(pwanted)-log(sommaperet));
        //calcolo i pesi finali usando l'et corretto
        for(int i=0;i<nwalkers;i++){
            pesi[i]*=exp(-tau/2.*(detemp[i]-2*et[countdmc]*npassidmc));
            outpesi[countdmc]+=pesi[i];
        }
        // birth-death
        std::vector<int> topop; //non posso pop-are mentre sono nel ciclo altrimenti perdo il conto
        int naggiunti = 0;
        for(int i=0;i<nwalkers;i++){
            if(pesi[i]>1){
                if((double)rand()/(double)RAND_MAX<(pesi[i]-1)){ //devo generare un nuovo walker
                    pesi[i]=1;
                    Walker * temp = new Walker(nparticelle, latoscatola, &parametripot[0], 2, &parametrilogprob[0], 2, walkers[i]->px, walkers[i]->py, walkers[i]->pz);
                    walkers.push_back(temp);
                    Observables * tempobs = new Observables(npassiblock);
                    obs.push_back(tempobs);
                    pesi.push_back(1);
                    naggiunti++;
                }
            }else{
                if((double)rand()/(double)RAND_MAX>pesi[i]){
                    topop.push_back(i);
                }
            }
        }
        nwalkers+=naggiunti;
        int nremove = 0;
        //pop-o quello che devo pop-are
        for(int i=0;i<topop.size();i++){
            // prima elimino l'oggetto salvato in quella parte del vector
            delete walkers[topop[i]];
            delete obs[topop[i]];
            for(int j=topop[i];j<nwalkers-1;j++){
                walkers[j] = walkers[j+1];
                obs[j] = obs[j+1];
                pesi[j] = pesi[j+1];
            }
            walkers.pop_back();
            obs.pop_back();
            pesi.pop_back();
            nremove++;
            // std::cout << "walker " << topop[i] << " removed at step " << countdmc << std::endl; 
            for(int j=i+1;j<topop.size();j++){
                topop[j]--;
            }
        }
        nwalkers -= nremove;

        std::cout << countdmc << " " << nwalkers << " " << naggiunti << " " << nremove << std::endl;
    }

    //stampo su file i valori di et: 
    FILE * fpet = fopen("mt_val/valet.dat", "w");
    for(int i=0;i<npassidmc;i++){
        fprintf(fpet, "%lf\n", et[i]);
    }
    fclose(fpet);
    //stampo su file i valori di et cumulata
    double * etcum = mediacumulata(et, npassidmc, 0);
    FILE * fpetcum=fopen("mt_val/valetcum.dat", "w");
    for(int i=0;i<npassidmc;i++){
        fprintf(fpetcum, "%lf\n", etcum[i]);
    }
    fclose(fpetcum);
    //stampo su file i valori di energia media sui walker
    FILE * fpen=fopen("mt_val/valen.dat", "w");
    for(int i=0;i<npassidmc*npassiblock;i++){
        fprintf(fpen, "%lf\n", energiemediebig[i]);
    }
    fclose(fpen);
    //stampo su file la media cumulata
    double * energiemediebigcum = mediacumulata(energiemediebig, npassidmc*npassiblock, 0);
    FILE * fpencum=fopen("mt_val/valencum.dat", "w");
    for(int i=0;i<npassidmc*npassiblock;i++){
        fprintf(fpencum, "%lf\n", energiemediebigcum[i]);
    }
    fclose(fpencum);

    FILE * fppesi=fopen("mt_val/outpesi.dat", "w");
    for(int i=0;i<npassidmc;i++){
        fprintf(fppesi, "%lf\n", outpesi[i]);
    }
    fclose(fppesi);
    
    FILE * fpken=fopen("mt_val/outken.dat", "w");
    for(int i=0;i<npassidmc;i++){
        fprintf(fpken, "%lf\n", kenout[i]);
    }
    fclose(fpken);
    
    FILE * fppot=fopen("mt_val/outpot.dat", "w");
    for(int i=0;i<npassidmc;i++){
        fprintf(fppot, "%lf\n", potout[i]);
    }
    fclose(fppot);
    
    FILE * fpstim1=fopen("mt_val/outstim1.dat", "w");
    for(int i=0;i<npassidmc;i++){
        fprintf(fpstim1, "%lf\n", stim1out[i]);
    }
    fclose(fpstim1);
    
    FILE * fpstim2=fopen("mt_val/outstim2.dat", "w");
    for(int i=0;i<npassidmc;i++){
        fprintf(fpstim2, "%lf\n", stim2out[i]);
    }
    fclose(fpstim2);
    //stampo a video media e varianza scema
    double * energiemediebigfinal = media(energiemediebig, npassidmc*npassiblock, 0);
    std::cout << "media e varianza energia " << energiemediebigfinal[0] << " " << energiemediebigfinal[1] << std::endl; 
    
    return 0;
}