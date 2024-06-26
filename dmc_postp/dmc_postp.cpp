#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>

using namespace std;

double * med(int nelem, double * vmed, double * vvar){
    double * out = new double[2];
    for(int i=0;i<nelem;i++){
        out[0] += vmed[i];
        out[1] += vvar[i];
    }
    out[0] /= nelem;
    out[1] /= nelem;
    return out;
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

int main(){
    string workdir="varmcnomrt_detbal";
    int nwalk=40;
    int npassi=500000;
    double tau=0.00005;
    double pwanted=nwalk;
    double varianzeken[nwalk]={1759.649067,1769.964963,1671.022573,1819.724379,1840.346912,1836.601824,1690.786235,1878.579833,1854.561585,1827.178519,1580.615617,1599.387445,1722.045567,1770.429627,1684.016113,1803.515849,1894.836537,1778.374706,1827.269007,1861.006164,1868.186680,1704.088783,1754.892440,1921.270107,1797.765767,1924.619162,1928.810656,1557.804859,1783.980356,1854.891621,834.362777,1703.506746,1840.563830,1965.581038,1042.439686,1748.537893,1635.848794,1868.521812,1612.406655,1828.258978};
    double varianzepot[nwalk]={68135.251261,71466.971890,67506.346697,69824.175301,73646.430680,72395.776634,69048.803611,72993.751431,72085.702149,74226.581905,64279.213983,66348.905244,68626.222035,70672.601696,66798.428089,68194.703546,72460.751074,71627.036161,71700.212310,71541.725366,74219.210511,67789.797998,69631.812826,74716.718691,71993.520066,78022.858218,76093.645115,65642.303383,67997.723167,72639.295934,30131.473666,67830.579476,74091.506727,79276.700898,37369.557899,69728.651984,64180.967719,75509.911467,65618.115291,70750.391881};

    ifstream * filesken;
    ifstream * filespot;
    ofstream outpesi("outpesi.dat");
    ofstream outvalori("outvalori.dat");
    filesken = new ifstream[nwalk];
    filespot = new ifstream[nwalk];
    for(int i=0;i<nwalk;i++){
        string nomefileken = "../"+workdir+"/mt_val_dim/ken_"+to_string(i)+".dat";
        string nomefilepot = "../"+workdir+"/mt_val_dim/pot_"+to_string(i)+".dat";
        filesken[i].open(nomefileken, std::ifstream::in);
        filespot[i].open(nomefilepot, std::ifstream::in);
    }
    int npassidmc = npassi-1;
    double * et = new double[npassidmc];
    double * etcum;
    double * medieE = new double[npassidmc];
    double * medieEcum;
    double ** medieEcumpesata;
    double * medieEtest = new double[npassidmc];
    double * varEtest = new double[npassidmc];
    double * varE = new double[npassidmc];
    double * kenpartenza = new double[nwalk];
    double * kenarrivo = new double[nwalk];
    double * potpartenza = new double[nwalk];
    double * potarrivo = new double[nwalk];
    double * deltae = new double[nwalk];
    double tempsommede=0;
    double * pesi = new double[nwalk];
    for(int i=0;i<nwalk;i++){
        pesi[i]=1;
    }

    for(int i=0;i<nwalk;i++){
        string kenstr;
        getline(filesken[i], kenstr);
        kenarrivo[i] = stod(kenstr);
        string potstr;
        getline(filespot[i], potstr);
        potarrivo[i] = stod(potstr);
    }

    for(int i=0;i<npassidmc;i++){
        tempsommede=0;
        for(int j=0;j<nwalk;j++){
            kenpartenza[j]=kenarrivo[j];
            potpartenza[j]=potarrivo[j];
            string kenstr;
            getline(filesken[j], kenstr);
            kenarrivo[j] = stod(kenstr);
            string potstr;
            getline(filespot[j], potstr);
            potarrivo[j] = stod(potstr);
            deltae[j] = potarrivo[j]+kenarrivo[j]+(potpartenza[j]+kenpartenza[j]);
            tempsommede+=pesi[j]*exp(-tau/2.*deltae[j]);
        }   

        et[i]=1./tau*(log(pwanted)-log(tempsommede));
        for(int j=0;j<nwalk;j++){
            pesi[j] *= exp(-tau/2.*(deltae[j]-2*et[i]));
            outpesi << pesi[j] << " ";
        }
        outpesi << endl;

        double kenmedtemp=0;
        double kenmedtemptest=0;
        double kenvartemptest=0;
        double potmedtemp=0;
        double potmedtemptest=0;
        double potvartemptest=0;
        double kenvartemp=0;
        double potvartemp=0;
        double pesopermedia=0;
        for(int j=0;j<nwalk;j++){
            kenmedtemp+=kenarrivo[j]*pesi[j];
            kenmedtemptest+=kenarrivo[j];
            potmedtemp+=potarrivo[j]*pesi[j];
            potmedtemptest+=potarrivo[j];
            kenvartemp+=varianzeken[j]*pesi[j];
            kenvartemptest+=varianzeken[j];
            potvartemp+=varianzepot[j]*pesi[j];
            potvartemptest+=varianzepot[j];
            pesopermedia+=pesi[j];
        }
        kenmedtemp/=pesopermedia;
        potmedtemp/=pesopermedia;
        kenvartemp/=pesopermedia;
        potvartemp/=pesopermedia;
        kenmedtemptest/=nwalk;
        potmedtemptest/=nwalk;
        kenvartemptest/=nwalk;
        potvartemptest/=nwalk;
        medieE[i] = kenmedtemp+potmedtemp;
        varE[i] = kenvartemp+potvartemp;
        medieEtest[i] = kenmedtemptest+potmedtemptest;
        varEtest[i] = kenvartemptest+potvartemptest;
        outvalori << et[i] << " " << medieE[i] << " " << varE[i] << endl;

        if(i%10000==0){
            cout << "\r " << i << flush;
        }
    }
    cout << endl;

    double * valfinale = med(npassidmc, medieE, varE);
    double * controllo = med(npassidmc, medieEtest, varEtest);

    int metapassidmc=(int)((double)npassidmc/2.+.5);

    double * valfinaleprimameta = med(metapassidmc, medieE, varE);
    double * valfinalesecondameta = med(metapassidmc, &medieE[metapassidmc], &varE[metapassidmc]);

    double etmedia;
    double etvar;
    medvar(et, npassidmc, &etmedia, &etvar);

    cout << "media e varianza con la media pesata" << valfinale[0] << " " << sqrt(valfinale[1]/npassidmc) << endl; 
    cout << "media e varianza con et " << etmedia << " " << sqrt(etvar/npassidmc) << endl; 
    cout << "controllo media e variaza " << controllo[0] << " " << sqrt(controllo[1]/npassidmc) << endl;
    cout << "medie e varianze prima meta: " << valfinaleprimameta[0] << " " << sqrt(valfinaleprimameta[1]/metapassidmc) << endl;
    cout << "medie e varianze seconda meta: " << valfinalesecondameta[0] << " " << sqrt(valfinalesecondameta[1]/metapassidmc) << endl;

    for(int i=0;i<nwalk;i++){
        filesken[i].close();
        filespot[i].close();
    }
    outpesi.close();

    etcum = mediacumulata(et, npassidmc, 0);
    medieEcum = mediacumulata(medieE, npassidmc, 0);
    ofstream outvaloricum("outvaloricum.dat");
    for(int i=0;i<npassidmc;i++){
        outvaloricum << etcum[i] << " " << medieEcum[i] << endl;
    }
    outvaloricum.close();


    delete [] et;
    delete [] medieE;
    delete [] medieEtest;
    delete [] varEtest;
    delete [] varE;
    delete [] kenpartenza;
    delete [] kenarrivo;
    delete [] potpartenza;
    delete [] potarrivo;
    delete [] deltae;
    delete [] pesi;
 
    return 0;
}