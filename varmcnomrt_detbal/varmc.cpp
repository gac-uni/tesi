#include "../libs/gpd.h"
#include "variables.h"

void invertivettori(double * v1, double * v2, int lunghezza){
    double temp;
    for(int i=0;i<lunghezza;i++){
        temp = v1[i];
        v1[i]=v2[i];
        v2[i]=temp;
    }
}

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
        int plotposizioni(std::string nomefile){
            GnuplotDriver posizionigp;
            posizionigp.ls("points");
            posizionigp.limx(-20,20);
            posizionigp.limy(-20,20);
            posizionigp.limz(-20,20);
            posizionigp.fpath(nomefile);
            posizionigp.fext("jpg");
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
        int nelem;
        Observables(int nelementi){
            nelem = nelementi;
            vPot = new double[nelem];
            vKen = new double[nelem];
            vStim1 = new double[nelem];
            vStim2 = new double[nelem];
            vEn = new double[nelem];
            vAcc = new double[nelem];
        }
        ~Observables(){
            delete [] vPot;
            delete [] vKen;
            delete [] vStim1;
            delete [] vStim2;
            delete [] vEn;
            delete [] vAcc;
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

int main(int argc, char ** argv){
    srand(69);

    switch(argc){
        case 1:
            break;
        case 2:
            std::sscanf(argv[1], "%lf", &tau);
            break;
        case 3:
            std::sscanf(argv[1], "%lf", &tau);
            std::sscanf(argv[2], "%d", &npassi);   
            npassiplot = npassi/freqcampionamento+0.5;         
            break;
    }
    
    GnuplotDriver kengp;
    GnuplotDriver potgp;
    GnuplotDriver kenditgp;
    GnuplotDriver potditgp;
    GnuplotDriver stimatore1gp;
    GnuplotDriver stimatore2gp;
    GnuplotDriver controllogp;
    GnuplotDriver energiagp;
    GnuplotDriver energiaditgp;
    GnuplotDriver pesoditgp;

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
    controllogp.fpath(plotdir+"posiniziale");
    controllogp.ls("points");
    controllogp.plot(posx, posy, posz, nparticelle);
    
    Walker walker(nparticelle, latoscatola, &parametripot[0], 2, &parametrilogprob[0], 2, posx, posy, posz);
    Observables obs(npassiplot);
    Observables obs_totale(npassi);

    // --- evoluzione temporale ---
    int contatoreplot=0;
    double naccettati=0;
    // vettori con le comp dei gradienti della configurazione di partenza che serve per calcolare la drift velocity
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
    for(int i=0; i<npassi; i++){
        // al passo iniziale devo calcolare manualmente i gradienti
        if(i==0){
            for(int j=0;j<nparticelle;j++){
                walker.computegradi(j, &gradlogprex[j], &gradlogprey[j], &gradlogprez[j]);
            }
        }
        // la funzione di montecarlo prende come input i gradienti dei punti in modo da non doverli ricalcolare
        int accettato;
        accettato=montecarlo(&walker, tau, gradlogprex, gradlogprey, gradlogprez, gradlogsuccx, gradlogsuccy, gradlogsuccz);
        naccettati+=accettato;
        // --- calcolo osservabili ---
        // la funzione dell'energia cinetica restituisce i vettori dei gradienti aggiornati alla nuova posizione, che al passo successivo saranno quelli della posizione di partenza
        double * tempken=walker.ken(gradlogsuccx, gradlogsuccy, gradlogsuccz);
        obs_totale.vKen[i]=tempken[0];
        obs_totale.vStim1[i]=tempken[1];
        obs_totale.vStim2[i]=tempken[2];
        obs_totale.vPot[i]=walker.pot();
        obs_totale.vEn[i]=tempken[0]+obs_totale.vPot[i];
        delete [] tempken;
        obs_totale.vAcc[i]=(double)naccettati/(double)(i+1);
        for(int j=0;j<nparticelle;j++){
            gradlogprex[j]=gradlogsuccx[j];
            gradlogprey[j]=gradlogsuccy[j];
            gradlogprez[j]=gradlogsuccz[j];
        }
        if(i%freqcampionamento==0){
            // --- se necessario salvo separatamente i dati per i plot rapidi ---
            obs.vPot[contatoreplot] = obs_totale.vPot[i];
            obs.vKen[contatoreplot] = obs_totale.vKen[i];
            obs.vStim1[contatoreplot] = obs_totale.vStim1[i];
            obs.vStim2[contatoreplot] = obs_totale.vStim2[i];
            obs.vEn[contatoreplot] = obs_totale.vEn[i];
            obs.vAcc[contatoreplot] = obs_totale.vAcc[i];
            // --- salvo le posizioni dei punti come plot ---
            // std::string nome = plotdir+"pos/pos";
            // nome = nome+std::to_string(contatoreplot);
            // walker.plotposizioni(nome);
            contatoreplot++;
        }
        if(i%(int)(npassi/100 + 0.5)==0){            
            std::cout<<"\rcompletamento: " << (double)i/(double)npassi*100. << "%, probabilità di accettazione: " << (double)naccettati/(double)i << std::flush;
        }
    }
    std::cout<<std::endl;

    // --- stampo su file tutte le osservabili complete --- 
    std::string nomeken=valdir+"ken.dat";
    std::string nomestim1=valdir+"stim1.dat";
    std::string nomestim2=valdir+"stim2.dat";
    std::string nomepot=valdir+"pot.dat";
    std::string nomeen=valdir+"en.dat";
    std::string nomeacc=valdir+"acc.dat";
    std::string nomepesi=valdir+"pesi.dat";
    FILE * fpken=fopen(nomeken.c_str(), "w");
    FILE * fpstim1=fopen(nomestim1.c_str(), "w");
    FILE * fpstim2=fopen(nomestim2.c_str(), "w");
    FILE * fppot=fopen(nomepot.c_str(), "w");
    FILE * fpen=fopen(nomeen.c_str(), "w");
    FILE * fpacc=fopen(nomeacc.c_str(), "w");
    for(int i=0;i<npassi;i++){
        fprintf(fpken, "%lf\n", obs_totale.vKen[i]);
        fprintf(fpstim1, "%lf\n", obs_totale.vStim1[i]);
        fprintf(fpstim2, "%lf\n", obs_totale.vStim2[i]);
        fprintf(fppot, "%lf\n", obs_totale.vPot[i]);
        fprintf(fpen, "%lf\n", obs_totale.vEn[i]);
        fprintf(fpacc, "%lf\n", obs_totale.vAcc[i]);
    }
    fclose(fpken);
    fclose(fpstim1);
    fclose(fpstim2);
    fclose(fppot);
    fclose(fpen);
    fclose(fpacc);

    // faccio le medie cumulate
    double * kencum;
    double * potcum;
    double * stimatore1cum;
    double * stimatore2cum;
    double * energiacum;
    kencum=mediacumulata(obs.vKen, npassiplot, 0);
    potcum=mediacumulata(obs.vPot, npassiplot, 0);
    stimatore1cum=mediacumulata(obs.vStim1, npassiplot, 0);
    stimatore2cum=mediacumulata(obs.vStim2, npassiplot, 0);
    energiacum=mediacumulata(obs.vEn, npassiplot, 0);

    // grafici
    // kengp.fpath(plotdir+"ken");
    // kengp.ls("line");
    // kengp.noprint();
    // kengp.plot(kencum, npassiplot, 1);

    // potgp.fpath(plotdir+"pot");
    // potgp.ls("line");
    // potgp.noprint();
    // potgp.plot(potcum, npassiplot, 1);

    kenditgp.fpath(plotdir+"kendit");
    kenditgp.ls("line");
    kenditgp.t("energia cinetica");
    kenditgp.y("energia [adimensionale]");
    kenditgp.plot(obs.vKen, npassiplot, 1);

    potditgp.fpath(plotdir+"potdit");
    potditgp.ls("line");
    potditgp.t("potenziale");
    potditgp.y("energia [adimensionale]");
    potditgp.plot(obs.vPot, npassiplot, 1);

    stimatore1gp.fpath(plotdir+"stimatore1");
    stimatore1gp.ls("line");
    stimatore1gp.plot(obs.vStim1, npassiplot, 1);

    stimatore2gp.fpath(plotdir+"stimatore2");
    stimatore2gp.ls("line");
    stimatore2gp.plot(obs.vStim2, npassiplot, 1);

    // energiagp.fpath(plotdir+"energia");
    // energiagp.ls("line");
    // energiagp.plot(energiacum, npassiplot, 1)

    energiaditgp.fpath(plotdir+"energiadit");
    energiaditgp.ls("line");
    energiaditgp.t("potenziale");
    energiaditgp.y("energia [adimensionale]");
    energiaditgp.plot(obs.vEn, npassiplot, 1);
    
    return 0;
}